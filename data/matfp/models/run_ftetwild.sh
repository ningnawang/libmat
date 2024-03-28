#!/bin/bash

# Set the number of threads to use
THREADS=4

# Find all files under the specified folder and run the command on each file in parallel
find . -type f -name "*.obj" | xargs -I {} -P $THREADS sh -c '/home/junanita/ninwang/fTetWild/build/FloatTetwild_bin -i "$1" -l 1' sh {}
