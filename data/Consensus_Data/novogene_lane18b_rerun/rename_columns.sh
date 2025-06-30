#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

input_file="$1"
output_file="$2"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Error: File '$input_file' not found."
  exit 1
fi

# Read the file and replace the column headers
awk 'BEGIN { IGNORECASE = 1 } \
    NR == 1 { \
        for (i = 1; i <= NF; i++) {
            if (tolower($i) ~ /ct_d[06]/) $i = "ct";
            if (tolower($i) ~ /depth_d[06]/) $i = "depth";
        }
    } \
    { print }' "$input_file" > "$output_file"

# Notify the user
echo "Processed file saved to '$output_file'."
