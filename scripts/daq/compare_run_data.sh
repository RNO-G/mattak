#!/bin/bash
#
# compare_station_files_with_summit.sh
#
# Description:
#   Compares a remote data directory with a local copy using rsync (no checksum check)
#   to identify identical and different directories.
#   Can optionally remove identical directories from the remote after confirmation.
#   Skips empty directories (on the remote), i.e., directories smaller than a
#   specified size limit (default: 24 KB).
#
#   Two modes:
#     default (station -> summit): remote s<station_id>:/data/daq
#                                  vs. local /data/ingress/station<station_id>/
#     --chicago-to-summit:         remote greenland:/data/archived/station<station_id>
#                                  vs. local /data/full/raw/station<station_id>/
#                                  (run this on the uchicago server)
#
# Usage:
#   $0 <station_id> [--size-limit <KB>] [--no-remove] [--chicago-to-summit] [--run-pattern <glob>]
#
# Arguments:
#   <station_id>         The station ID
#   --size-limit <KB>    Minimum directory size to include in comparison (default: 24 KB)
#   --no-remove          Skip automatic removal of identical directories (default: will remove after confirmation)
#   --chicago-to-summit  Compare summit archive (greenland:/data/archived/station<id>)
#                        against local uchicago copy (/data/full/raw/station<id>).
#                        Always a dry run: implies --no-remove.
#   --run-pattern <glob> Only consider remote directories matching this glob (default: *)
#
# Workflow:
#   1. Connects to the remote host and lists directories in the remote dir (skipping those smaller than size limit)
#   3. Compares each remote directory with its local counterpart using rsync
#   4. Categorizes directories as IDENTICAL or DIFFERENT or MISSING LOCALLY
#   5. Displays summary with counts and total sizes
#   6. If --no-remove is not set, prompts user to remove identical directories from remote
#
# Example:
#   ./compare_station_files_with_summit.sh 13
#   ./compare_station_files_with_summit.sh 13 --size-limit 5120 --no-remove
#   ./compare_station_files_with_summit.sh 13 --chicago-to-summit --run-pattern 'run*'
#

STATION_ID="$1"
if [[ -z "$STATION_ID" ]]; then
    echo "Usage: $0 <station_id> [--size-limit <KB>] [--no-remove] [--chicago-to-summit] [--run-pattern <glob>]"
    echo "  --size-limit: Size limit in KB (default: 24)"
    echo "  --no-remove: Skip removal of identical directories on remote (needs confirmation anyway...)"
    echo "  --chicago-to-summit: Compare greenland:/data/archived/station<id> with local /data/full/raw/station<id> (implies --no-remove)"
    echo "  --run-pattern: Only consider remote directories matching this glob (default: *)"
    exit 1
fi

SIZE_LIMIT_KB=24  # Default size limit in KB
REMOVE_FLAG=true
CHICAGO_MODE=false
RUN_PATTERN="*"

# Parse additional arguments
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        --size-limit) SIZE_LIMIT_KB="$2"; shift ;;
        --no-remove) REMOVE_FLAG=false ;;
        --chicago-to-summit) CHICAGO_MODE=true ;;
        --run-pattern) RUN_PATTERN="$2"; shift ;;
        *) ;;
    esac
    shift
done

if [[ "$CHICAGO_MODE" == true ]]; then
    LOCAL_DIR="/data/full/raw/station${STATION_ID}"
    REMOTE_DIR="/data/archived/station${STATION_ID}"
    HOST="greenland"
    REMOVE_FLAG=false  # chicago mode is always a dry run: never remove from the summit archive
else
    LOCAL_DIR="/data/ingress/station${STATION_ID}"
    REMOTE_DIR="/data/daq"
    HOST="s${STATION_ID}"
fi

# Verify we can connect to the host before doing any work
# (BatchMode prevents hanging on a password prompt)
echo "Checking connection to $HOST..."
if ! ssh -q -o BatchMode=yes -o ConnectTimeout=10 "$HOST" true; then
    echo "Error: Cannot connect to host '$HOST' (check ssh config, network, or keys)"
    exit 1
fi

# Verify local directory exists
if [[ ! -d "$LOCAL_DIR" ]]; then
    echo "Error: Local directory $LOCAL_DIR does not exist"
    exit 1
fi

skip_run=""
if [[ "$CHICAGO_MODE" == false ]]; then
    # Query the current run from the station and subtract 1 to get the run to skip
    echo "Querying current run from $HOST..."
    next_run=$(ssh -q "$HOST" "cat /rno-g/var/runfile" 2>/dev/null | tr -d '[:space:]')
    if [[ -z "$next_run" ]]; then
        echo "Warning: Could not read /rno-g/var/runfile from $HOST"
    else
        skip_run=$((next_run - 1))
        echo "Current run: $skip_run, will skip ..."
    fi
fi

# Get remote directory sizes into associative array
echo "Connecting to $HOST and analyzing $REMOTE_DIR..."
declare -A REMOTE_DIR_SIZES
declare -a REMOTE_DIRS  # to have sorted list of directories with size > limit
while IFS=$'\t' read -r size dir; do
    if (( size > SIZE_LIMIT_KB )); then
        REMOTE_DIR_SIZES["$dir"]=$size
        REMOTE_DIRS+=("$dir")
    fi
done < <(ssh -q "$HOST" "cd $REMOTE_DIR && du -sk $RUN_PATTERN 2>/dev/null" | sort -k2)

if [[ ${#REMOTE_DIRS[@]} -eq 0 ]]; then
    echo "No remote directories matching '$RUN_PATTERN' with size > ${SIZE_LIMIT_KB} KB found in $HOST:$REMOTE_DIR"
    exit 0
fi

# Compare directories
echo ""
echo "Comparison of directories with size > ${SIZE_LIMIT_KB} KB:"
echo "=========================================="

# Initialize arrays for results
identical_dirs=()
different_dirs=()
identical_size=0
different_size=0

# Create temporary file with list of directories to sync (only those above size limit)
rsync_includes=$(mktemp)
trap "rm -f '$rsync_includes'" EXIT
printf '%s\n' "${REMOTE_DIRS[@]}" > "$rsync_includes"

# Use single rsync call to compare only directories above size limit (much faster)
echo "Running rsync comparison on directories with size > ${SIZE_LIMIT_KB} KB..."
rsync_output=$(rsync -anr --itemize-changes --files-from="$rsync_includes" -e "ssh -q" "${HOST}:${REMOTE_DIR}/" "${LOCAL_DIR}/" 2>&1)
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

# Parse rsync output to identify which directories have differences
declare -A dir_has_differences
echo "Comparing rsync output..."
# Extract directories with changes from rsync output
while IFS= read -r line; do
    # Skip empty lines and metadata lines (sent, received, total)
    if [[ -z "$line" || "$line" =~ ^(sent|received|total) ]]; then
        continue
    fi

    # rsync format: "changeinfo path" - extract the path part
    if [[ $line =~ ^[^[:space:]]+[[:space:]]+(.+)$ ]]; then
        filepath="${BASH_REMATCH[1]}"
        # Get top-level directory (run folder)
        top_dir=$(echo "$filepath" | cut -d'/' -f1)
        [[ -n "$top_dir" ]] && dir_has_differences["$top_dir"]=true
    fi
done <<< "$rsync_output"

# Categorize directories based on rsync findings
for dir in "${REMOTE_DIRS[@]}"; do
    # Skip the run specified by skip_run if it matches the directory name
    if [[ -n "$skip_run" ]]; then
        dir_num="${dir#run}"  # Extract numeric part by removing "run" prefix
        if [[ "$dir_num" == "$skip_run" ]]; then
            echo "Skip current run: $dir"
            continue
        fi
    fi

    dir_size=${REMOTE_DIR_SIZES["$dir"]}
    local_path="${LOCAL_DIR}/${dir}/"

    # Check if local directory exists
    if [[ ! -d "$local_path" ]]; then
        echo "MISSING LOCALLY: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
    elif [[ "${dir_has_differences[$dir]}" == "true" ]]; then
        echo "DIFFERENT: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
    else
        echo "IDENTICAL: $dir"
        identical_dirs+=("$dir")
        ((identical_size += dir_size))
    fi
done

# Extract network traffic from rsync's own summary line ("sent X bytes  received Y bytes ...")
traffic_line=$(awk '/^sent /{gsub(/,/,""); printf "%.2f MB sent + %.2f MB received = %.2f MB", $2/1048576, $5/1048576, ($2+$5)/1048576}' <<< "$rsync_output")

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
    echo ""
    echo "The following command will be executed to remove directories from $HOST:$REMOTE_DIR:"
    echo "Command: ssh $HOST rm -rf $(printf "${REMOTE_DIR}/%s " "${identical_dirs[@]}")"
    read -p "Are you sure you want to remove these directories? (yes/no): " confirm

    if [[ "$confirm" == "yes" ]]; then
        echo "Removing identical directories from remote..."

        # Execute single SSH command with all directories at once
        ssh -q "$HOST" "rm -rf $(printf "${REMOTE_DIR}/%s " "${identical_dirs[@]}")" 2>&1
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
fi
