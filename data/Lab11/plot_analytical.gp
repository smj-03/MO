set title "Rozwiazanie Analityczne"
set xlabel "x"
set ylabel "U(x, t)"
set grid
set key outside
plot  "data/Lab11/analytical_values.txt" using 1:2 with lines, \
      "data/Lab11/explicit_values.txt" using 1:2 with points