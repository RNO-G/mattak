# Ported from NuRadioReco/utilities/logging.py

import logging

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
    level : int, default=25
        The logging level to use for the base logger

    """
    logger = logging.getLogger(name)

    if len(logger.handlers) > 0:  # method hasHandlers() also checks parents -> ends up at root logger
        # Don't change the logger if it already exists
        logger.warning(f"Logger {name} already has handlers. Not changing anything, returning the existing logger...")
        return logger
    logger.propagate = False

    # Create a StreamHandler with fancy formatter
    handler = logging.StreamHandler()
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

            self.FORMATS = {
                logging.DEBUG: grey + "%(levelname)s - " + reset + format,
                logging.INFO: green + "%(levelname)s - " + reset + format,
                logging.WARNING: purple + "%(levelname)s - " + reset + format,
                logging.ERROR: red + "%(levelname)s - " + reset + format,
                logging.CRITICAL: red + "%(levelname)s - " + reset + format
            }

        def format(self, record):
            log_fmt = self.FORMATS.get(record.levelno)
            formatter = logging.Formatter(log_fmt)
            return formatter.format(record)


    formatter = CustomFormatter(
        format='\033[33m%(asctime)s - \033[32m%(name)s - \033[0m%(message)s',
        datefmt="%H:%M:%S"
    )

    return formatter

_setup_logger()