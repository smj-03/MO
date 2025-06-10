set xlabel "log10(h)"
set ylabel "log10(|error|)"
set grid
set key top right

plot 'data/Lab11/errors_explicit.txt' using (log10($1)):(log10($2)) with lines title 'KMB error', \
 'data/Lab11/errors_implicit_thomas.txt' using (log10($1)):(log10($2)) with lines title 'THOMAS error', \
  'data/Lab11/errors_implicit_lu.txt' using (log10($1)):(log10($2)) with lines title 'LU error'