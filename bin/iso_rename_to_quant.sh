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

find "$directory" -type f -name "*.sf.gz" -exec rename 's/\.sf\.gz$/quant.sf.gz/' {} +