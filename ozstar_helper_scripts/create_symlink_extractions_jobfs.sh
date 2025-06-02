#!/bin/bash

# Check if the field name argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <fieldname>"
    exit 1
fi

# Assign input parameter to a variable
FIELDNAME="$1"

# Define source and destination directories using the input parameter
SOURCE_DIR="/fred/oz041/lkawin/NISPureParallel/FIELDS/${FIELDNAME}/Extractions"
#DEST_DIR="/fred/oz041/lkawin/NISPureParallel/FIELDS/${FIELDNAME}/RedshiftFitting"
DEST_DIR="$JOBFS/data/lkawin/FIELDS/${FIELDNAME}/RedshiftFitting"
# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Navigate to the destination directory
cd "$DEST_DIR" || { echo "Failed to change directory to $DEST_DIR"; exit 1; }

# Create symbolic links for specific file patterns
find "$SOURCE_DIR" -type f \( \
    -name "*full.fits" -o \
    -name "*-extracted.fits" -o \
    -name "*1D.fits" -o \
    -name "*row.fits" -o \
    -name "*stack.pkl" \
\) -exec ln -s {} . \;

echo "Symbolic links created successfully in $DEST_DIR"

