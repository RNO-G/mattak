#!/bin/bash
#
# compare_run_data.sh
#
# Description:
#   Compares a source data directory with a copy of it using rsync (no checksum check)
#   to identify identical and different directories. The source is either on a remote
#   host (default / --chicago-to-summit) or on the same host (--summit).
#   Can optionally remove identical directories from the source after confirmation.
#   Skips empty directories (on the source), i.e., directories smaller than a
#   specified size limit (default: 24 KB).
#   Only directories named "run*" are ever compared and hence ever removed.
#
#   Three modes:
#     default (station -> summit): remote s<station_id>:/data/daq
#                                  vs. local /data/ingress/station<station_id>/
#     --chicago-to-summit:         remote greenland:/data/archived/station<station_id>
#                                  vs. local /data/full/raw/station<station_id>/
#                                  (run this on the uchicago server)
#     --summit (ingress -> archive): local /data/ingress/station<station_id>
#                                  vs. local /data/archived/station<station_id>
#                                  (run this on the summit server, no ssh involved)
#
# Usage:
#   compare_run_data.sh <station_id> [--size-limit <KB>] [--no-remove] [--auto-approve-remove] [--chicago-to-summit] [--summit] [--host <host>] [--run-pattern <glob>]
#
# Arguments:
#   <station_id>         The station ID
#   -h, --help           Show usage and exit
#   --size-limit <KB>    Minimum directory size to include in comparison (default: 24 KB)
#   --no-remove          Skip automatic removal of identical directories (default: will remove after confirmation)
#   --auto-approve-remove  Remove identical directories without asking for confirmation.
#                        Mutually exclusive with --no-remove (and hence --chicago-to-summit).
#   --chicago-to-summit  Compare summit archive (greenland:/data/archived/station<id>)
#                        against local uchicago copy (/data/full/raw/station<id>).
#                        Always a dry run: implies --no-remove.
#   --summit             Compare local /data/ingress/station<id> against local
#                        /data/archived/station<id> (both on the summit server).
#   --host <host>        Host holding the source data, overriding the per-mode default
#                        (s<station_id>, or "greenland" with --chicago-to-summit).
#                        With --summit the source is local unless --host is given.
#   --run-pattern <glob> Only consider source directories matching this glob (default: *).
#                        Entries not named "run*" are ignored regardless of the pattern.
#   --verbose            Print the raw rsync output (for debugging)
#
# Workflow:
#   1. Lists the run directories in the source dir (skipping non-"run*" entries and
#      those smaller than size limit), connecting to the remote host first if the
#      source is remote
#   2. Asks the station (s<station_id>) for the run currently being taken and excludes it
#      (default and --summit mode; best effort, a station that is offline is not an error)
#   3. Compares each source directory with its counterpart in the copy using rsync
#   4. Categorizes directories as IDENTICAL or DIFFERENT or MISSING IN COPY
#   5. Displays summary with counts and total sizes
#   6. If --no-remove is not set, prompts user to remove identical directories from the source
#      (no prompt with --auto-approve-remove)
#   7. In the dry-run case (--no-remove / --chicago-to-summit), writes the identical runs
#      (station<id>/run<id>, one per line) to a text file for later manual removal
#
# Example:
#   ./compare_run_data.sh 13
#   ./compare_run_data.sh 13 --size-limit 5120 --no-remove
#   ./compare_run_data.sh 13 --auto-approve-remove
#   ./compare_run_data.sh 13 --chicago-to-summit --run-pattern 'run*'
#   ./compare_run_data.sh 13 --summit
#   ./compare_run_data.sh 13 --host s13-lte
#

usage() {
    echo "Usage: $0 <station_id> [--size-limit <KB>] [--no-remove] [--auto-approve-remove] [--chicago-to-summit] [--summit] [--host <host>] [--run-pattern <glob>]"
    echo "  -h, --help: Show this help message and exit"
    echo "  --size-limit: Size limit in KB (default: 24)"
    echo "  --no-remove: Skip removal of identical directories in the source (needs confirmation anyway...)"
    echo "  --auto-approve-remove: Remove identical directories without confirmation (mutually exclusive with --no-remove)"
    echo "  --chicago-to-summit: Compare greenland:/data/archived/station<id> with local /data/full/raw/station<id> (implies --no-remove)"
    echo "  --summit: Compare local /data/ingress/station<id> with local /data/archived/station<id> (no ssh)"
    echo "  --host: Host holding the source data (default: s<station_id>, or 'greenland' with --chicago-to-summit)"
    echo "  --run-pattern: Only consider source directories matching this glob (default: *, non-'run*' entries are always ignored)"
    echo "  --verbose: Print the raw rsync output (for debugging)"
}

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

STATION_ID="$1"
if [[ -z "$STATION_ID" ]]; then
    usage
    exit 1
fi

SIZE_LIMIT_KB=24  # Default size limit in KB
REMOVE_FLAG=true
AUTO_APPROVE_REMOVE=false
CHICAGO_MODE=false
SUMMIT_MODE=false
RUN_PATTERN="*"
VERBOSE=false
HOST_OVERRIDE=""

# Parse additional arguments
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        --size-limit) SIZE_LIMIT_KB="$2"; shift ;;
        --no-remove) REMOVE_FLAG=false ;;
        --auto-approve-remove) AUTO_APPROVE_REMOVE=true ;;
        --chicago-to-summit) CHICAGO_MODE=true ;;
        --summit) SUMMIT_MODE=true ;;
        --host) HOST_OVERRIDE="$2"; shift ;;
        --run-pattern) RUN_PATTERN="$2"; shift ;;
        --verbose) VERBOSE=true ;;
        -h|--help) usage; exit 0 ;;
        *) ;;
    esac
    shift
done

if [[ "$CHICAGO_MODE" == true && "$SUMMIT_MODE" == true ]]; then
    echo "Error: --chicago-to-summit and --summit are mutually exclusive"
    exit 1
fi

# --chicago-to-summit implies --no-remove, so it conflicts with --auto-approve-remove as well
if [[ "$AUTO_APPROVE_REMOVE" == true && ( "$REMOVE_FLAG" == false || "$CHICAGO_MODE" == true ) ]]; then
    echo "Error: --auto-approve-remove and --no-remove (implied by --chicago-to-summit) are mutually exclusive"
    exit 1
fi

if [[ "$AUTO_APPROVE_REMOVE" == true ]]; then
    echo "################################################################################"
    echo "#                                  WARNING                                     #"
    echo "#  --auto-approve-remove is set: identical run directories will be REMOVED     #"
    echo "#  from the source WITHOUT asking for confirmation. Press Ctrl-C to abort.     #"
    echo "################################################################################"
fi

# SOURCE_DIR/LOCAL_DIR denote the source of the data and the copy we compare it
# against. In --summit mode both live on this host, hence SOURCE_IS_REMOTE.
# LOCAL_LABEL identifies the host holding the copy (used for the dry-run filename).
SOURCE_IS_REMOTE=true
if [[ "$CHICAGO_MODE" == true ]]; then
    LOCAL_DIR="/data/full/raw/station${STATION_ID}"
    SOURCE_DIR="/data/archived/station${STATION_ID}"
    HOST="greenland"
    LOCAL_LABEL="uchicago"
    REMOVE_FLAG=false  # chicago mode is always a dry run: never remove from the summit archive
elif [[ "$SUMMIT_MODE" == true ]]; then
    LOCAL_DIR="/data/archived/station${STATION_ID}"
    SOURCE_DIR="/data/ingress/station${STATION_ID}"
    HOST="s${STATION_ID}"  # the data is local, the host is only used to query the current run
    LOCAL_LABEL="summit_archived"
    SOURCE_IS_REMOTE=false
else
    LOCAL_DIR="/data/ingress/station${STATION_ID}"
    SOURCE_DIR="/data/daq"
    HOST="s${STATION_ID}"
    LOCAL_LABEL="summit"
fi

# An explicitly given host replaces the per-mode default. Passing it in --summit
# mode means the source dir is read over ssh from that host instead of locally.
if [[ -n "$HOST_OVERRIDE" ]]; then
    HOST="$HOST_OVERRIDE"
fi

if [[ "$SOURCE_IS_REMOTE" == true ]]; then
    SOURCE_DESC="${HOST}:${SOURCE_DIR}"
else
    SOURCE_DESC="$SOURCE_DIR"
fi

# True for plain "run*" directory names. Used as a safe guard on both ends of the
# script: nothing else is compared, and nothing else can be handed to rm -rf.
is_run_dir() {
    [[ "$1" == run* && "$1" != */* ]]
}

# Run a shell command on the side holding the source data: over ssh if that side
# is a remote host, in a local subshell otherwise (--summit).
source_exec() {
    if [[ "$SOURCE_IS_REMOTE" == true ]]; then
        ssh -q "$HOST" "$1"
    else
        bash -c "$1"
    fi
}

if [[ "$SOURCE_IS_REMOTE" == true ]]; then
    # Verify we can connect to the host before doing any work
    # (BatchMode prevents hanging on a password prompt)
    echo "Checking connection to $HOST..."
    if ! ssh -q -o BatchMode=yes -o ConnectTimeout=10 "$HOST" true; then
        echo "Error: Cannot connect to host '$HOST' (check ssh config, network, or keys)"
        exit 1
    fi
elif [[ ! -d "$SOURCE_DIR" ]]; then
    echo "Error: Source directory $SOURCE_DIR does not exist"
    exit 1
fi

# Verify the directory holding the copy exists
if [[ ! -d "$LOCAL_DIR" ]]; then
    echo "Error: Local directory $LOCAL_DIR does not exist"
    exit 1
fi

skip_run=""
if [[ "$CHICAGO_MODE" == false ]]; then
    # Query the next (/rno-g/var/runfile) run from the station and subtract 1 to get the run to skip.
    # The query is best effort: a station that is currently unreachable must not abort the
    # comparison (in --summit mode we have not talked to it before), hence BatchMode/ConnectTimeout
    # so we neither hang on a password prompt nor on a dead link.
    echo "Querying current run from $HOST..."
    next_run=$(ssh -q -o BatchMode=yes -o ConnectTimeout=30 "$HOST" "cat /rno-g/var/runfile" 2>/dev/null | tr -d '[:space:]')
    if [[ -z "$next_run" ]]; then
        echo "Warning: Could not read /rno-g/var/runfile from $HOST"
    else
        skip_run=$((next_run - 1))
        echo "Current run: $skip_run, will skip ..."
    fi
fi

# Get source directory sizes into associative array
echo "Analyzing $SOURCE_DESC..."
declare -A SOURCE_DIR_SIZES
declare -a SOURCE_DIRS  # to have sorted list of directories with size > limit
while IFS=$'\t' read -r size dir; do
    # Safe guard: only run directories are ever compared (and hence ever removed).
    # This drops stray files/directories in the source (e.g. logs, tarballs) that
    # --run-pattern would otherwise let through. is_run_dir also rejects names
    # containing a "/" so nothing outside SOURCE_DIR can end up in the rm command.
    if ! is_run_dir "$dir"; then
        [[ "$VERBOSE" == true ]] && echo "Skipping non-run entry: $dir"
        continue
    fi

    if (( size > SIZE_LIMIT_KB )); then
        SOURCE_DIR_SIZES["$dir"]=$size
        SOURCE_DIRS+=("$dir")
    fi
done < <(source_exec "cd $SOURCE_DIR && du -sk $RUN_PATTERN 2>/dev/null" | sort -k2)

if [[ ${#SOURCE_DIRS[@]} -eq 0 ]]; then
    echo "No directories matching '$RUN_PATTERN' with size > ${SIZE_LIMIT_KB} KB found in $SOURCE_DESC"
    exit 0
fi

echo "Found ${#SOURCE_DIRS[@]} directories matching '$RUN_PATTERN' with size > ${SIZE_LIMIT_KB} KB in $SOURCE_DESC"

# Initialize arrays for results
identical_dirs=()
different_dirs=()
identical_size=0
different_size=0

# Create temporary file with list of directories to sync (only those above size limit)
rsync_includes=$(mktemp)
trap "rm -f '$rsync_includes'" EXIT
printf '%s\n' "${SOURCE_DIRS[@]}" > "$rsync_includes"

# Use single rsync call to compare only directories above size limit (much faster)
# --exclude drops rsync's own leftover temp files from failed prior transfers,
# named "<original-name>.XXXXXX" with a leading dot and a 6-char random suffix
# (e.g. .000172.wf.dat.gz.KsUBEc), so they don't cause spurious DIFFERENT results.
echo "Running rsync comparison on ${#SOURCE_DIRS[@]} directories with size > ${SIZE_LIMIT_KB} KB..."
if [[ "$SOURCE_IS_REMOTE" == true ]]; then
    rsync_source="${HOST}:${SOURCE_DIR}/"
    rsync_transport=(-e "ssh -q")
else
    # Local source: no transport, rsync compares the two paths directly
    rsync_source="${SOURCE_DIR}/"
    rsync_transport=()
fi
rsync_output=$(rsync -anr --itemize-changes --stats --exclude='.*.??????' --files-from="$rsync_includes" "${rsync_transport[@]}" "$rsync_source" "${LOCAL_DIR}/" 2>&1)
rsync_exit_code=$?

# Check if rsync failed
if [[ $rsync_exit_code -ne 0 ]]; then
    echo "Error: rsync command failed with exit code $rsync_exit_code. Aborting."
    exit 1
fi

# Check if rsync output is empty
if [[ -z "$rsync_output" ]]; then
    echo "Error: rsync output is empty. Aborting."
    exit 1
fi

if [[ "$VERBOSE" == true ]]; then
    echo "---- raw rsync output ----"
    echo "$rsync_output"
    echo "---------------------------"
fi

# Parse rsync output to identify which directories have differences
declare -A dir_has_differences
echo "Comparing rsync output..."
# Extract directories with changes from rsync output
while IFS= read -r line; do
    # Skip empty lines, the --stats block, and the final summary (sent, received, total)
    if [[ -z "$line" || "$line" =~ ^(sent|received|total|Number|Total|Literal|Matched|Unmatched|File) ]]; then
        continue
    fi

    # rsync format: "itemcode path" - extract the itemcode and path
    if [[ $line =~ ^([^[:space:]]+)[[:space:]]+(.+)$ ]]; then
        itemcode="${BASH_REMATCH[1]}"
        filepath="${BASH_REMATCH[2]}"
        # itemcode[1] is the file type (f=file, d=dir, L=symlink, ...). Directory
        # entries only ever reflect attribute noise (e.g. mtime) since their
        # content is captured by the individual file entries within them, so
        # they don't indicate an actual data difference and would otherwise
        # cause every run to be falsely marked DIFFERENT.
        [[ "${itemcode:1:1}" == "d" ]] && continue
        # Get top-level directory (run folder)
        top_dir=$(echo "$filepath" | cut -d'/' -f1)
        [[ -n "$top_dir" ]] && dir_has_differences["$top_dir"]=true
    fi
done <<< "$rsync_output"

# Categorize directories based on rsync findings
for dir in "${SOURCE_DIRS[@]}"; do
    # Skip the run specified by skip_run if it matches the directory name
    if [[ -n "$skip_run" ]]; then
        dir_num="${dir#run}"  # Extract numeric part by removing "run" prefix
        if [[ "$dir_num" == "$skip_run" ]]; then
            [[ "$VERBOSE" == true ]] && echo "Skip current run: $dir"
            continue
        fi
    fi

    dir_size=${SOURCE_DIR_SIZES["$dir"]}
    local_path="${LOCAL_DIR}/${dir}/"

    # Check if local directory exists
    if [[ ! -d "$local_path" ]]; then
        [[ "$VERBOSE" == true ]] && echo "MISSING IN COPY: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
    elif [[ "${dir_has_differences[$dir]}" == "true" ]]; then
        [[ "$VERBOSE" == true ]] && echo "DIFFERENT: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
    else
        [[ "$VERBOSE" == true ]] && echo "IDENTICAL: $dir"
        identical_dirs+=("$dir")
        ((identical_size += dir_size))
    fi
done

# Extract network traffic from rsync's own summary line ("sent X bytes  received Y bytes ...").
# Only meaningful when the source is remote; a local comparison never hits the network.
traffic_line=""
if [[ "$SOURCE_IS_REMOTE" == true ]]; then
    traffic_line=$(awk '/^sent /{gsub(/,/,""); printf "%.2f MB sent + %.2f MB received = %.2f MB", $2/1048576, $5/1048576, ($2+$5)/1048576}' <<< "$rsync_output")
fi

# Print summary
echo ""
echo "=========================================="
echo "SUMMARY:"
echo "  Identical: ${#identical_dirs[@]} dirs ($(echo "scale=2; $identical_size / 1048576" | bc) GB) - "
echo "$(IFS=', '; echo "${identical_dirs[*]}")"
echo "  Different: ${#different_dirs[@]} dirs ($(echo "scale=2; $different_size / 1048576" | bc) GB) - "
echo "$(IFS=', '; echo "${different_dirs[*]}")"
[[ -n "$traffic_line" ]] && echo "  Network traffic (rsync): $traffic_line"
echo "=========================================="

# Handle removal of identical directories if requested
if [[ "$REMOVE_FLAG" == true && ${#identical_dirs[@]} -gt 0 ]]; then
    # Safe guard (second line of defence, the listing is filtered already): never
    # build an rm command from anything that is not a plain run directory
    for dir in "${identical_dirs[@]}"; do
        if ! is_run_dir "$dir"; then
            echo "Error: refusing to remove '$dir': not a run directory. Aborting."
            exit 1
        fi
    done

    rm_cmd="rm -rf $(printf "${SOURCE_DIR}/%s " "${identical_dirs[@]}")"
    echo ""
    echo "The following command will be executed to remove directories from $SOURCE_DESC:"
    if [[ "$SOURCE_IS_REMOTE" == true ]]; then
        echo "Command: ssh $HOST $rm_cmd"
    else
        echo "Command: $rm_cmd"
    fi
    if [[ "$AUTO_APPROVE_REMOVE" == true ]]; then
        echo "--auto-approve-remove is set: removing without confirmation!"
        confirm="yes"
    else
        read -p "Are you sure you want to remove these directories? (yes/no): " confirm
    fi

    if [[ "$confirm" == "yes" ]]; then
        echo "Removing identical directories from $SOURCE_DESC..."

        # Execute a single command with all directories at once
        source_exec "$rm_cmd" 2>&1
        if [[ $? -eq 0 ]]; then
            echo "  ✓ Successfully removed all identical directories"
        else
            echo "  ✗ Error during removal"
        fi
        echo "Removal complete."
    else
        echo "Removal cancelled."
    fi
elif [[ "$REMOVE_FLAG" == true ]]; then
    echo "No identical directories to remove."
elif [[ ${#identical_dirs[@]} -gt 0 ]]; then
    # Dry run: write the identical runs to a file for later manual removal.
    # Metadata (host holding the copy, station, date) goes into the filename
    # so the file itself contains only paths, one per line.
    RUNLIST_FILE="identical_runs_${LOCAL_LABEL}_station${STATION_ID}_$(date +%Y-%m-%d).txt"
    printf "station${STATION_ID}/%s\n" "${identical_dirs[@]}" > "$RUNLIST_FILE"

    echo ""
    echo "Dry run: wrote ${#identical_dirs[@]} identical runs to $RUNLIST_FILE"
    echo "To remove them, run on the host holding the data (from the directory containing station${STATION_ID}/):"
    echo "  for d in \$(cat $RUNLIST_FILE); do rm -vr \"\$d\"; done"
fi
