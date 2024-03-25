gnuplot << EOF
set terminal png size 1920,1080
set output "Q.png"

set xlabel "Time (seconds)"
set ylabel "Q (flow)"

plot "Q_plots_s.txt" with lines title 'Start' linecolor rgb "red", "Q_plots_m.txt" with lines title 'Middle' linecolor rgb "blue", "Q_plots_e.txt" with lines title 'End' linecolor rgb "green"
EOF

gnuplot << EOF
set terminal png size 1920,1080
set output "P.png"

set xlabel "Time (seconds)"
set ylabel "P (pressure)"

plot "P_plots_s.txt" with lines title 'Start' linecolor rgb "red", "P_plots_m.txt" with lines title 'Middle' linecolor rgb "blue", "P_plots_e.txt" with lines title 'End' linecolor rgb "green"
EOF