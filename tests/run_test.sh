#!/bin/bash
path="$HOME/Luka_work/code_area/Open-Reading-Frame-Finder/results/benchmark_results"
if [ ! -d "$path" ]; then
    cd ../results/
    mkdir benchmark_results
    echo "Created results folder"
else
    echo "Results folder exists, output will be redirected"
fi

./benchmark_test_script.py > "$path/benchmark_output.txt"
