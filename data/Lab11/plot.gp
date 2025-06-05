set term qt 0
set title "Rozwiazanie Analityczne"
set xlabel "x"
set ylabel "U(x, t)"
set grid
set key outside
plot for [i=2:7] "data/Lab11/wykres_analityczny" using 1:i with lines title sprintf("t=%.1f", (i-2)*0.1)

set term qt 1
set title "Rozwiazanie Numeryczne"
set xlabel "x"
set ylabel "U(x, t)"
set grid
set key outside
plot for [i=2:7] "data/Lab11/wykres_numeryczny" using 1:i with lines title sprintf("t=%.1f", (i-2)*0.1)