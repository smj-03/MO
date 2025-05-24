set multiplot layout 1,3 title "Wykresy funkcji x/(1+x^4) oraz jej interpolacji metoda Lagrange'a"

set xlabel "x"
set ylabel "f(x)"
set grid
plot  "data/Lab12/punkty_funkcji.txt" using 1:2 with lines linewidth 3 linecolor rgb "green" title "x/(1+x^4)", \
      "data/Lab12/punkty_interpolacji_lagrange.txt" using 1:2 with lines linecolor rgb "red" title "Interpolacja (baza Lagrange'a)", \
      "data/Lab12/punkty_interpolacji_czebyszew.txt" using 1:2 with lines linecolor rgb "blue" title "Interpolacja (baza Czebyszew'a)"

set xlabel "x"
set ylabel "f(x)"
set grid
plot  "data/Lab12/punkty_funkcji.txt" using 1:2 with lines linewidth 3 linecolor rgb "green" title "x/(1+x^4)", \
      "data/Lab12/punkty_interpolacji_lagrange.txt" using 1:2 with lines linecolor rgb "red" title "Interpolacja (baza Lagrange'a)"

set xlabel "x"
set ylabel "f(x)"
set grid
plot  "data/Lab12/punkty_funkcji.txt" using 1:2 with lines linewidth 3 linecolor rgb "green" title "x/(1+x^4)", \
      "data/Lab12/punkty_interpolacji_czebyszew.txt" using 1:2 with lines linecolor rgb "blue" title "Interpolacja (baza Czebyszew'a)"

unset multiplot