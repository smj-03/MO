set title "Wykres funkcji x/(1+x^4)"
set xlabel "x"
set ylabel "f(x)"
set grid
set term wxt
plot  "data/Lab12/punkty_funkcji.txt" using 1:2 with lines linewidth 3 title "Wykres funkcji x/(1+x^4)", \