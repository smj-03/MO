set title "Wykres funkcji x/(1+x^4)"
set grid
set term wxt
plot \
"data/Lab2/wyniki_1.txt" using 1:2 with lines title "log10(x) od log10(error)", \
"data/Lab2/wyniki_2.txt" using 1:2 with lines title "log10(x) od log10(alt_error)", \
-16 with lines lc rgb "red" title "Blad reprezentacji"