#!/bin/bash

# config
MAX_MEM_GB=8
MAX_TIME_SEC=30
ALGORITHM="-astar" # {-bfs | -idfs | -astar | -idastar | -gbfs}
INPUT="instances/15puzzle_instances.txt"
OUTPUT="results/astar15.csv"

MAX_MEM_KB=$((MAX_MEM_GB * 1024 * 1024))
> "$OUTPUT"

lineno=0

while IFS= read -r instance; do
    lineno=$((lineno + 1))
    echo "==> $lineno: $instance | $ALGORITHM"
    output=$( 
        ulimit -v "$MAX_MEM_KB"
        ulimit -t "$MAX_TIME_SEC"
        exec "./main" "$ALGORITHM" $instance
    ) 2>&1
    exit_status=$?
    if [ $exit_status -eq 0 ]; then
        echo "$lineno,$output" >> "$OUTPUT"
    else
        echo "$lineno,-,-,-,-,-" >> "$OUTPUT"
    fi
done < "$INPUT"
