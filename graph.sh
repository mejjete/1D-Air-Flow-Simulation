Q_s="${1}_Qplots_s.txt"
Q_m="${1}_Qplots_m.txt"
Q_e="${1}_Qplots_e.txt"

P_s="${1}_Pplots_s.txt"
P_m="${1}_Pplots_m.txt"
P_e="${1}_Pplots_e.txt"

gnuplot << EOF
set terminal png size 1920,1080
set output "Q.png"

set xlabel "Time (seconds)"
set ylabel "Q (flow)"

plot "$Q_s" with lines title 'Start' linecolor rgb "red", "$Q_m" with lines title 'Middle' linecolor rgb "blue", "$Q_e" with lines title 'End' linecolor rgb "green"
EOF

gnuplot << EOF
set terminal png size 1920,1080
set output "P.png"

set xlabel "Time (seconds)"
set ylabel "P (pressure)"

plot "$P_s" with lines title 'Start' linecolor rgb "red", "$P_m" with lines title 'Middle' linecolor rgb "blue", "$P_e" with lines title 'End' linecolor rgb "green"
EOF
