#!/bin/bash

# Define the field range, excluding sex-07
for i in {11..32}; do
    if [ "$i" == "22" ]; then
        continue  # Skip sex-22
    fi
    field="sex-$i"

    FIELD_DIR="/fred/oz041/lkawin/NISPureParallel/FIELDS/$field"
    BACKUP_DIR="${FIELD_DIR}_backup"
    TAR_FILE="$FIELD_DIR/${field}_RedshiftFitting.tar.gz"
    REDSHIFT_DIR="$FIELD_DIR/RedshiftFitting"
    TEMP_UNPACK_DIR="$FIELD_DIR/temp_unpack"
    SQUASHFS_FILE="${field}.sqfs"

    echo "Processing $field..."

    # Remove existing RedshiftFitting directory if it exists
    if [ -d "$REDSHIFT_DIR" ]; then
        echo "Removing existing $REDSHIFT_DIR..."
        rm -rf "$REDSHIFT_DIR"
    fi

    # Create RedshiftFitting directory
    mkdir -p "$REDSHIFT_DIR"

    # Unpack the tar.gz file if it exists
    if [ -f "$TAR_FILE" ]; then
        echo "Extracting $TAR_FILE..."

        # Create temporary unpack directory
        mkdir -p "$TEMP_UNPACK_DIR"

        # Extract contents to the temporary directory
        tar -xzf "$TAR_FILE" -C "$TEMP_UNPACK_DIR"

        # Check if extraction resulted in an extra nested directory structure
        if [ -d "$TEMP_UNPACK_DIR/$field/RedshiftFitting" ]; then
            echo "Detected nested directory structure, correcting it..."
            mv "$TEMP_UNPACK_DIR/$field/RedshiftFitting"/* "$REDSHIFT_DIR"
        else
            echo "Moving extracted files directly..."
            mv "$TEMP_UNPACK_DIR"/* "$REDSHIFT_DIR"
        fi

        # Remove temporary unpack directory
        rm -rf "$TEMP_UNPACK_DIR"

    else
        echo "Warning: $TAR_FILE not found. Skipping extraction."
    fi

    # Create SquashFS file for the entire field directory
    echo "Creating SquashFS for $field..."
    mksquashfs "$FIELD_DIR" "$SQUASHFS_FILE"

    # Rename FIELD_DIR to FIELD_DIR_backup after successful squashfs creation
    echo "Renaming $FIELD_DIR to $BACKUP_DIR..."
    mv "$FIELD_DIR" "$BACKUP_DIR"

    # Optionally remove the original FIELD_DIR (commented out for now)
    # echo "Removing $FIELD_DIR..."
    # rm -rf "$FIELD_DIR"

    echo "Finished processing $field."
done
