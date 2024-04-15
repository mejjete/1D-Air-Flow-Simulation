#!/bin/bash

arg="$1"
number=""
# Check if the argument starts with "E" and has more than one character
if [[ $arg =~ ^V([0-9]+)* ]]; then
    number="${BASH_REMATCH[1]}"
else
    echo "Invalid format: $arg"
    exit 1
fi

V_plots="${1}_Pplots.txt"
V_out="${1}_P.png"

gnuplot << EOF
set terminal png size 1920,1080
set output "$V_out"

set xlabel "Time (seconds)"
set ylabel "V_P$number (vertex pressure)"

plot "$V_plots" with lines title 'Start' linecolor rgb "red"
EOF
