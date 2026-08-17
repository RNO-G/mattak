# Ported from NuRadioReco/utilities/logging.py

import logging
from collections import Counter, defaultdict
from collections.abc import Mapping


class DeduplicatingStreamHandler(logging.StreamHandler):
    """
    A `logging.StreamHandler` which emits each distinct message only once.

    Records at or above `dedup_level` are identified by their *unformatted* message
    template (together with the logger name and level). The first record of each such
    group is emitted normally; subsequent ones are counted silently and their arguments
    are stored. A summary of all suppressed records is emitted by `emit_summary`, which
    is called automatically when the handler is closed (i.e. by `logging.shutdown` at
    interpreter exit).

    Records below `dedup_level`, and records carrying ``extra={"no_dedup": True}``, are
    passed through unchanged.

    Notes
    -----
    Because messages are keyed on the unformatted template, deferred formatting is
    required for deduplication to take effect: ``logger.warning("ch %d saturated", ch)``
    collapses into a single group, while ``logger.warning(f"ch {ch} saturated")``
    produces a distinct group per channel.

    Logging with *named* fields, i.e. a single dict as the only argument
    (``logger.warning("ev %(event)d bad", {"event": ev, "run": run})``), gives the summary
    field names to work with: it then groups the suppressed records by `group_keys` (e.g.
    per run) and reports the value range of the remaining fields per group. Fields do not
    have to appear in the message; the dict may carry extra context purely for the summary.

    The arguments of all suppressed records are stored, but the summary only reports the
    `max_groups` most affected groups and collapses field values to their range. Enable
    `logging.DEBUG` on the respective logger to list every group and value individually.

    Parameters
    ----------
    stream : file-like, optional
        Stream to write to. Passed to `logging.StreamHandler`; defaults to `sys.stderr`.
    dedup_level : int, default=logging.WARNING
        Minimum level at which records are deduplicated. Records below this level are
        always emitted.
    group_keys : sequence of str, default=("station", "run")
        Fields of dict-style arguments by which the summary sub-groups the suppressed
        records. Missing fields are ignored.
    max_groups : int, default=5
        Maximum number of such sub-groups reported individually (the most affected ones).
        The remaining ones are collapsed into a single count.
    """

    def __init__(self, stream=None, dedup_level=logging.WARNING, group_keys=("station", "run"),
                 max_groups=5):
        super().__init__(stream)
        self.dedup_level = dedup_level
        self.group_keys = group_keys
        self.max_groups = max_groups
        self._counts = Counter()
        self._samples = defaultdict(list)
        self._summarized = False

    def emit(self, record):
        """
        Emit a record, or count it as a repeat of one already emitted.

        Parameters
        ----------
        record : logging.LogRecord
            The record to handle.
        """
        if record.levelno < self.dedup_level or getattr(record, "no_dedup", False):
            super().emit(record)
            return

        key = (record.name, record.levelno, record.msg)
        if self._counts[key] == 0:
            super().emit(record)
        elif record.args:
            self._samples[key].append(record.args)
        self._counts[key] += 1

    def emit_summary(self):
        """
        Emit one summary line per group with suppressed records, then reset the state.

        Groups are reported in order of decreasing occurrence count; groups seen only
        once are skipped, having already been emitted in full. Calling this at a natural
        batch boundary (per run, per file) reports the summary while it is still in
        context and starts the next batch with clean counters. This is useful in
        long-running processes and interactive sessions, where the handler is never
        closed and the automatic summary would otherwise never appear.

        Field values are collapsed to their range, unless the group's logger is enabled for
        `logging.DEBUG`, in which case they are listed individually.
        """
        for key, n in sorted(self._counts.items(), key=lambda kv: -kv[1]):
            if n <= 1:
                continue
            name, levelno, msg = key
            samples = self._samples[key]
            verbose = logging.getLogger(name).isEnabledFor(logging.DEBUG)

            detail = self._format_samples(samples, verbose)
            if detail:
                detail = "\nThe suppressed calls had the following arguments:" + detail
                if not verbose and len(samples) > 1:
                    detail += "\n\t  (enable debug logging to list all arguments)"

            super().emit(logging.LogRecord(
                name=name, level=levelno, pathname=__file__, lineno=0,
                msg="\nThe following message was suppressed %d further time(s):\n\t%r%s",
                args=(n - 1, msg, detail), exc_info=None,
            ))

        self._counts.clear()
        self._samples.clear()

    def _format_samples(self, samples, verbose):
        """
        Render the stored arguments of one message group as indented summary lines.

        Dict-style arguments are sub-grouped by `group_keys` and reported for the
        `max_groups` most affected groups; positional arguments are listed as they are.

        Parameters
        ----------
        samples : list
            The stored `logging.LogRecord.args` of the suppressed records: either tuples
            (positional) or dicts (named).
        verbose : bool
            If True, list all values individually instead of collapsing them to a range.

        Returns
        -------
        detail : str
            Possibly empty (records without arguments).
        """
        if not samples:
            return ""

        # Positional args offer no field names to group by, so just list them (unwrapping
        # single-element tuples, i.e. "1, 2" rather than "(1,), (2,)").
        if not all(isinstance(sample, Mapping) for sample in samples):
            return "\n\t  args: " + self._format_values(
                [s[0] if len(s) == 1 else s for s in samples], verbose)

        groups = defaultdict(list)
        for sample in samples:
            groups[tuple((k, sample[k]) for k in self.group_keys if k in sample)].append(sample)

        # Most affected groups first. The sort is stable, so groups with an equal number of
        # occurrences stay in the order they were encountered (i.e. usually run order).
        ranked = sorted(groups.items(), key=lambda kv: -len(kv[1]))
        shown = ranked if verbose else ranked[:self.max_groups]

        detail = ""
        for group, entries in shown:
            fields = {k: [e[k] for e in entries if k in e]
                      for k in dict.fromkeys(k for e in entries for k in e) if k not in dict(group)}
            label = " | ".join(f"{k}={v}" for k, v in group)
            stats = " | ".join(f"{k}: {self._format_values(v, verbose)}" for k, v in fields.items())
            detail += f"\n\t  {label + ' ' if label else ''}({len(entries)}x){' ' + stats if stats else ''}"

        dropped = ranked[len(shown):]
        if dropped:
            keys = "/".join(k for k, _ in dropped[0][0]) or "group"
            detail += (f"\n\t  ... and {len(dropped)} further {keys} "
                       f"({sum(len(entries) for _, entries in dropped)}x in total)")

        return detail

    def _format_values(self, values, verbose):
        """
        Collapse the values of a single field to their range, e.g. the span of event numbers
        covered by the suppressed records. Non-numeric fields are reported by their number
        of distinct values instead.

        Parameters
        ----------
        values : list
        verbose : bool
            If True, list all values individually instead of collapsing them.

        Returns
        -------
        formatted : str
        """
        if verbose:
            return ", ".join(str(v) for v in values)

        if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values):
            return str(values[0]) if min(values) == max(values) else f"{min(values)}..{max(values)}"

        # str() so that unhashable values (numpy arrays, ...) can be counted as well
        unique = set(str(v) for v in values)
        return unique.pop() if len(unique) == 1 else f"{len(unique)} distinct values"

    def close(self):
        """
        Emit the summary (once) and close the handler.

        Called by `logging.shutdown` at interpreter exit. Exceptions raised while
        building the summary are swallowed, so that a failure here cannot disrupt
        interpreter shutdown.
        """
        if not self._summarized:
            self._summarized = True
            try:
                self.emit_summary()
            except Exception:
                pass
        super().close()

def _setup_logger(name="mattak"):
    """
    Set up the parent logger which all module loggers should pass their logs on to. If this one already
    exists, nothing is done and the logger is returned as is. Otherwise, a single new `logging.StreamHandler()`
    with a custom formatter is added.

    Notes
    -----
    This function is only meant to be called once, on import, as part of the `__init__.py`.

    Parameters
    ----------
    name : str, default="mattak"
        The name of the base logger

    """
    logger = logging.getLogger(name)

    if len(logger.handlers) > 0:  # method hasHandlers() also checks parents -> ends up at root logger
        # Don't change the logger if it already exists
        logger.warning(f"Logger {name} already has handlers. Not changing anything, returning the existing logger...")
        return logger
    logger.propagate = False

    # Create a StreamHandler with fancy formatter
    handler = DeduplicatingStreamHandler()
    handler.setFormatter(get_fancy_formatter())
    handler.setLevel(1)  # we want the handler to be accepting all records from child loggers

    # Then add our custom handler to the logger
    logger.addHandler(handler)

    return logger


def get_fancy_formatter():
    """
    Returns the formatter used in the NuRadio logger.

    Returns
    -------
    formatter : logging.Formatter
    """

    class CustomFormatter(logging.Formatter):

        def __init__(self, format, datefmt):
            super().__init__(datefmt=datefmt)
            grey = "\033[38;1m"
            yellow = "\033[33;1m"
            purple = "\033[35;1m"
            green = "\033[32;1m"
            red = "\033[31;1m"
            reset = "\033[0m"

            # One formatter per level, built once: constructing them in format() would
            # drop `datefmt` and allocate a Formatter per record.
            self.FORMATS = {
                level: logging.Formatter(color + "%(levelname)s - " + reset + format, datefmt=datefmt)
                for level, color in [
                    (logging.DEBUG, grey),
                    (logging.INFO, green),
                    (logging.WARNING, purple),
                    (logging.ERROR, red),
                    (logging.CRITICAL, red),
                ]
            }
            # Used for custom levels (e.g. 25), which have no color assigned
            self.default_format = logging.Formatter(format, datefmt=datefmt)

        def format(self, record):
            return self.FORMATS.get(record.levelno, self.default_format).format(record)


    formatter = CustomFormatter(
        format='\033[33m%(asctime)s - \033[32m%(name)s - \033[0m%(message)s',
        datefmt="%H:%M:%S"
    )

    return formatter

_setup_logger()