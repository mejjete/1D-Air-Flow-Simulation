#!/bin/bash

arg="$1"
number=""
# Check if the argument starts with "E" and has more than one character
if [[ $arg =~ ^E([0-9]+)* ]]; then
    number="${BASH_REMATCH[1]}"
else
    echo "Invalid format: $arg"
    exit 1
fi

Q_s="${1}_Qplots_s.txt"
Q_m="${1}_Qplots_m.txt"
Q_e="${1}_Qplots_e.txt"
Q_out="${1}_Q.png"

P_s="${1}_Pplots_s.txt"
P_m="${1}_Pplots_m.txt"
P_e="${1}_Pplots_e.txt"
P_out="${1}_P.png"

if ! test -f "$Q_s"; then
    echo "$Q_s does not exist."
    exit 1
fi

if ! test -f "$Q_m"; then
    echo "$Q_m does not exist."
    exit 1
fi

if ! test -f "$Q_e"; then
    echo "$Q_e does not exist."
    exit 1
fi


if ! test -f "$P_s"; then
    echo "$P_s does not exist."
    exit 1
fi

if ! test -f "$P_m"; then
    echo "$P_m does not exist." 
    exit 1
fi

if ! test -f "$P_e"; then
    echo "$P_e does not exist."
    exit 1
fi

gnuplot << EOF
set terminal png size 1920,1080
set output "$Q_out"

set xlabel "Time (seconds)"
set ylabel "Q$number (flow)"

plot "$Q_s" with lines title 'Start' linecolor rgb "red", "$Q_m" with lines title 'Middle' linecolor rgb "blue", "$Q_e" with lines title 'End' linecolor rgb "green"
EOF

gnuplot << EOF
set terminal png size 1920,1080
set output "$P_out"

set xlabel "Time (seconds)"
set ylabel "P$number (pressure)"

plot "$P_s" with lines title 'Start' linecolor rgb "red", "$P_m" with lines title 'Middle' linecolor rgb "blue", "$P_e" with lines title 'End' linecolor rgb "green"
EOF
