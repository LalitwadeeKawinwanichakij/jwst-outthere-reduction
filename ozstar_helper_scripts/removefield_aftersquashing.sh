#!/bin/bash

FILELIST="tobesquash.dat"

if [ ! -f "$FILELIST" ]; then
  echo "Error: $FILELIST not found!"
  exit 1
fi

while read -r FIELDNAME; do
  [[ -z "$FIELDNAME" || "$FIELDNAME" == \#* ]] && continue

  SQUASH_FILE="/fred/oz408/NISPureParallel/squashfs_FIELDS/${FIELDNAME}.sqsh"
  FIELD_DIR="/fred/oz408/NISPureParallel/FIELDS/${FIELDNAME}"

  if [ -f "$SQUASH_FILE" ]; then
    echo "✔ Found squash file: $SQUASH_FILE"
    if [ -d "$FIELD_DIR" ]; then
      echo "⚠ Directory exists: $FIELD_DIR"
      echo -n "Do you want to remove the entire directory $FIELD_DIR? This will NOT prompt for every file. [y/N] " </dev/tty
      read -r confirm </dev/tty
      if [[ "$confirm" =~ ^[Yy]$ ]]; then
        rm -r "$FIELD_DIR"
        echo "✅ Removed: $FIELD_DIR"
      else
        echo "⏭ Skipped: $FIELD_DIR"
      fi
    else
      echo "ℹ Directory does not exist: $FIELD_DIR"
    fi
  else
    echo "✘ Squash file not found: $SQUASH_FILE"
  fi

  echo "-----------------------------"

done < "$FILELIST"
