#!/bin/bash

# This script compares the contents of two directories, checking for missing files and size differences.
# dir2 can have additional files, but all files in dir1 must be present in dir2 with the same size.
# For example you can check whether all files in INGRESS are also present in ARCHIVED

# Usage: ./check_dir1_in_dir2.sh dir1 dir2

# Remove trailing slashes if any
DIR1="${1%/}"
DIR2="${2%/}"

if [ -z "$DIR1" ] || [ -z "$DIR2" ]; then
    echo "Usage: $0 <dir1> <dir2>"
    exit 1
fi

# Recursively find all files in DIR1
find "$DIR1" -type f | while read -r file1; do
    # relative path of file in dir1
    rel="${file1#$DIR1/}"
    file2="$DIR2/$rel"

    if [ ! -f "$file2" ]; then
        echo "Missing in DIR2: $rel"
    else
        if [ "$(stat -c%s "$file1")" -ne "$(stat -c%s "$file2")" ]; then
            echo "Files differ (size): $rel"
	fi
    fi
done

echo "Check done."
