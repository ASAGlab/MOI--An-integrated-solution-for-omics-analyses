#!/bin/bash

# Check if the directory path argument is provided
if [[ -z "$1" ]]; then
  echo "Please provide the directory path as an argument."
  exit 1
fi

directory="$1"

# Check if the directory exists
if [[ ! -d "$directory" ]]; then
  echo "Directory '$directory' does not exist."
  exit 1
fi

# Change to the directory
cd "$directory"

# Iterate over each file in the directory
for file in *; do
  # Check if the item is a file
  if [[ -f "$file" ]]; then
    # Extract the filename without the extension
    filename="${file%.*}"
    
    # Remove the .sf.gz extension from the filename
    filename="${file%.sf.gz}"
    
    # Create a subdirectory with the same name as the file (if it doesn't exist)
    mkdir -p "$filename"
    
    # Move the file to the subdirectory and rename
    mv "$file" "$filename/quant.sf.gz"
    # unzip
    gunzip "$filename/quant.sf.gz"

  fi
done
