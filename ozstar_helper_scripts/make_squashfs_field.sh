#!/bin/bash

# Path to the field list
#FIELD_LIST="FIELDS/fields_dor_to_hya.txt"
FIELD_LIST="tobesquash.dat"
# Base directories
BASE_DIR="/fred/oz408/NISPureParallel/FIELDS"
LOG_DIR="/fred/oz408/NISPureParallel/logs"
SQSH_DIR="/fred/oz408/NISPureParallel/squashfs_FIELDS"
# Loop through each field in the list
while read -r field; do
    # Skip empty lines
    [[ -z "$field" ]] && continue

    #LOG_FILE="${LOG_DIR}/${field}_pipeline.log"

    # Check if log file exists
    #if [[ ! -f "$LOG_FILE" ]]; then
    #    echo "Log file not found for $field, skipping."
    #    continue
    #fi

    # Check if last line of log file contains "Pipeline complete ..."
    #if tail -n 1 "$LOG_FILE" | grep -q "Pipeline complete ..."; then
    echo "Creating SquashFS for $field..."
    FIELD_DIR="${BASE_DIR}/${field}"
    SQUASHFS_FILE="${SQSH_DIR}/${field}.sqsh"

    mksquashfs "$FIELD_DIR" "$SQUASHFS_FILE"
    #else
    #    echo "Pipeline not complete for $field, skipping."
    #fi

done < "$FIELD_LIST"

