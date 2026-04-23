#!/bin/bash
#
# compare_station_files_with_summit.sh
#
# Description:
#   Compares stations' /data/daq directory content with the local /data/ingress/station<station_id>/
#   Uses rsync without checksum check to identify identical and different directories.
#   Can optionally remove identical directories from the remote station after confirmation.
#   Skips empty directories (on the remote station), i.e., directories smaller than a
#   specified size limit (default: 24 KB).
#
# Usage:
#   $0 <station_id> [--size-limit <KB>] [--no-remove]
#
# Arguments:
#   <station_id>         The station ID (will connect to s<station_id>)
#   --size-limit <KB>    Minimum directory size to include in comparison (default: 24 KB)
#   --no-remove          Skip automatic removal of identical directories (default: will remove after confirmation)
#
# Workflow:
#   1. Connects to remote host s<station_id> and lists directories in /data/daq (skipping those smaller than size limit)
#   3. Compares each remote directory with its local counterpart under /data/ingress/station<station_id>/ using rsync
#   4. Categorizes directories as IDENTICAL or DIFFERENT or MISSING LOCALLY
#   5. Displays summary with counts and total sizes
#   6. If --no-remove is not set, prompts user to remove identical directories from remote
#
# Example:
#   ./compare_station_files_with_summit.sh 13
#   ./compare_station_files_with_summit.sh 13 --size-limit 5120 --no-remove
#

STATION_ID="$1"
if [[ -z "$STATION_ID" ]]; then
    echo "Usage: $0 <station_id> [--size-limit <KB>] [--no-remove]"
    echo "  --size-limit: Size limit in KB (default: 24)"
    echo "  --no-remove: Skip removal of identical directories on remote (needs confirmation anyway...)"
    exit 1
fi

SIZE_LIMIT_KB=24  # Default size limit in KB
REMOVE_FLAG=true

# Parse additional arguments
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        --size-limit) SIZE_LIMIT_KB="$2"; shift ;;
        --no-remove) REMOVE_FLAG=false ;;
        *) ;;
    esac
    shift
done

LOCAL_DIR="/data/ingress/station${STATION_ID}"

# Verify local directory exists
if [[ ! -d "$LOCAL_DIR" ]]; then
    echo "Error: Local directory $LOCAL_DIR does not exist"
    exit 1
fi

HOST="s${STATION_ID}"

# Get remote directory sizes into associative array
echo "Connecting to $HOST and analyzing /data/daq..."
declare -A REMOTE_DIR_SIZES
declare -a REMOTE_DIRS  # to have sorted list of directories with size > limit
while IFS=$'\t' read -r size dir; do
    if (( size > SIZE_LIMIT_KB )); then
        REMOTE_DIR_SIZES["$dir"]=$size
        REMOTE_DIRS+=("$dir")
    fi
done < <(ssh -q "$HOST" "cd /data/daq && du -sk * 2>/dev/null" | sort -k2)

# Compare directories
echo ""
echo "Comparison of directories with size > ${SIZE_LIMIT_KB} KB:"
echo "=========================================="

# Initialize arrays for results
identical_dirs=()
different_dirs=()
identical_size=0
different_size=0

# Compare each remote directory with local using rsync
for dir in "${REMOTE_DIRS[@]}"; do
    remote_path="${HOST}:/data/daq/${dir}/"
    local_path="${LOCAL_DIR}/${dir}/"
    dir_size=${REMOTE_DIR_SIZES["$dir"]}

    # Check if local directory exists
    if [[ ! -d "$local_path" ]]; then
        echo "MISSING LOCALLY: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
        continue
    fi

    # Use rsync with dry-run to compare
    # rsync -van --itemize-changes "$remote_path" "$local_path"
    diff_output=$(rsync -an --itemize-changes -e "ssh -q" "$remote_path" "$local_path" 2>&1)
    if [[ -z "$diff_output" ]]; then
        echo "IDENTICAL: $dir"
        identical_dirs+=("$dir")
        ((identical_size += dir_size))
    else
        echo "DIFFERENT: $dir"
        different_dirs+=("$dir")
        ((different_size += dir_size))
    fi
done

# Print summary
echo ""
echo "=========================================="
echo "SUMMARY:"
echo "  Identical: ${#identical_dirs[@]} dirs ($(echo "scale=2; $identical_size / 1048576" | bc) GB) - "
echo "$(IFS=', '; echo "${identical_dirs[*]}")"
echo "  Different: ${#different_dirs[@]} dirs ($(echo "scale=2; $different_size / 1048576" | bc) GB) - "
echo "$(IFS=', '; echo "${different_dirs[*]}")"
echo "=========================================="

# Handle removal of identical directories if requested
if [[ "$REMOVE_FLAG" == true && ${#identical_dirs[@]} -gt 0 ]]; then
    echo ""
    echo "The following directories will be removed from $HOST:/data/daq:"
    echo "Command:ssh $HOST rm -rf $(printf '/data/daq/%s ' "${identical_dirs[@]}")"
    read -p "Are you sure you want to remove these directories? (yes/no): " confirm

    if [[ "$confirm" == "yes" ]]; then
        echo "Removing identical directories from remote..."

        # Execute single SSH command with all directories at once
        ssh "$HOST" "rm -rf $(printf '/data/daq/%s ' "${identical_dirs[@]}")" 2>&1
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
