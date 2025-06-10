set title "Rozwiazanie Analityczne"
set xlabel "x"
set ylabel "U(x, t)"
set grid
set key outside
plot  "data/Lab11/analytical_values.txt" using 1:2 with lines, \
      "data/Lab11/explicit_values.txt" using 1:2 with points, \
      "data/Lab11/implicit_thomas_values.txt" using 1:2 with points pt 7 pointsize 0.3, \
      "data/Lab11/implicit_lu_values.txt" using 1:2 with points pt 4 pointsize 0.3, \