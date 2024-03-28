#!/bin/bash

# List of filenames to delete
files_to_delete=(
  "0000250_partstudio_02_model_ste_00_2048"
  "0002470_partstudio_04_model_ste_00_2048"
  "0004091_partstudio_00_model_ste_00_2048"
  "0004381_partstudio_07_model_ste_00_2048"
  "0006750_partstudio_00_model_ste_00_2048"
  "0011037_partstudio_01_model_ste_00_2048"
  "0011379_partstudio_03_model_ste_00_2048"
  "0011379_partstudio_24_model_ste_00_2048"
  "0011805_partstudio_00_model_ste_00_2048"
  "0011925_partstudio_00_model_ste_00_2048"
  "0014046_partstudio_01_model_ste_00_2048"
  "0015820_partstudio_01_model_ste_00_2048"
  "0000964_partstudio_02_model_ste_00_2048"
  "0007233_partstudio_00_model_ste_00_2048"
  "0009413_partstudio_00_model_ste_00_2048"
  "0012216_partstudio_00_model_ste_00_2048"
  "0012733_partstudio_00_model_ste_00_2048"
  "0005675_partstudio_03_model_ste_00_2048"
  "0011002_partstudio_00_model_ste_00_2048"
  "0013922_partstudio_02_model_ste_00_2048"
  "0001174_partstudio_01_model_ste_00_2048"
  "0002728_partstudio_01_model_ste_00_2048"
  "0010376_partstudio_11_model_ste_00_2048"
  "0004427_partstudio_00_model_ste_00_2048"
  "0010972_partstudio_00_model_ste_00_2048"
  "0015581_partstudio_01_model_ste_00_2048"
)

# Iterate over the filenames and delete the corresponding files with any extension
for filename in "${files_to_delete[@]}"; do
  # Use the 'find' command to locate and delete files
  find . -type f -name "${filename}.*" -exec rm {} \;
  echo "Deleted: ${filename}.*"
done

