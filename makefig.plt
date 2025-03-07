set term pdfcairo enhanced lw 1  ps 0.5 font ',12' size 7,7
set output "Num_test.pdf"
set encoding utf8
unset key

set multiplot layout 4,4

# First row: rho
set title "{\U+0394}x=0.02 {\U+03BE}=2.5 CFL=0.1"
set xlabel "x"
set ylabel "{\U+03C1}"
set yrange [-0.1:1.1]
plot 'out_sod_dense.txt' u 1:2 w l lw 2 lc "red" title 'Dense', \
     'out_sod_exact.txt' u 1:2 w l lw 0.5 lc -1 title 'Exact'

set xlabel "x"
set title "{\U+0394}x=0.12 {\U+03BE}=1 CFL=0.5"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:2 w p pt 6 lc "blue" title 'Sparse Large', \
     'out_sod_exact.txt' u 1:2 w l lw 0.5 lc -1 title 'Exact'

set xlabel "x"
set title "{\U+0394}x=0.02 {\U+03BE}=1 CFL=0.25"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:5 w p pt 6 lc "green" title 'Sparse Middle', \
     'out_sod_exact.txt' u 1:2 w l lw 0.5 lc -1 title 'Exact'

set xlabel "x"
set title "{\U+0394}x=0.02 {\U+03BE}=1 CFL=0.1"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:8 w p pt 6 lc "purple" title 'Sparse Low', \
     'out_sod_exact.txt' u 1:2 w l lw 0.5 lc -1 title 'Exact'

unset title

# Second row: v
set ylabel "v"
set yrange [-0.1:1.1]
plot 'out_sod_dense.txt' u 1:3 w l lw 2 lc "red" title 'Dense', \
     'out_sod_exact.txt' u 1:3 w l lw 0.5 lc -1 title 'Exact'

set ylabel "v"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:3 w p pt 6 lc "blue" title 'Sparse Large', \
     'out_sod_exact.txt' u 1:3 w l lw 0.5 lc -1 title 'Exact'

set ylabel "v"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:6 w p pt 6 lc "green" title 'Sparse Middle', \
     'out_sod_exact.txt' u 1:3 w l lw 0.5 lc -1 title 'Exact'

set ylabel "v"
set yrange [-0.1:1.1]
plot 'out_sod_sparse.txt' u 1:9 w p pt 6 lc "purple" title 'Sparse Low', \
     'out_sod_exact.txt' u 1:3 w l lw 0.5 lc -1 title 'Exact'

# Third row: p
set ylabel "p"
plot 'out_sod_dense.txt' u 1:4 w l lw 2 lc "red" title 'Dense', \
     'out_sod_exact.txt' u 1:4 w l lw 0.5 lc -1 title 'Exact'

set ylabel "p"
plot 'out_sod_sparse.txt' u 1:4 w p pt 6 lc "blue" title 'Sparse Large', \
     'out_sod_exact.txt' u 1:4 w l lw 0.5 lc -1 title 'Exact'

set ylabel "p"
plot 'out_sod_sparse.txt' u 1:7 w p pt 6 lc "green" title 'Sparse Middle', \
     'out_sod_exact.txt' u 1:4 w l lw 0.5 lc -1 title 'Exact'

set ylabel "p"
plot 'out_sod_sparse.txt' u 1:10 w p pt 6 lc "purple" title 'Sparse Low', \
     'out_sod_exact.txt' u 1:4 w l lw 0.5 lc -1 title 'Exact'

# Fourth row: s
set ylabel "s"
set yrange [-0.1:1.5]
s_f(p,r) = p/(r**1.4)-1

plot 'out_sod_dense.txt' u 1:(s_f($4,$2)) w l lw 2 lc "red" title 'Dense', \
     'out_sod_exact.txt' u 1:(s_f($4,$2)) w l lw 0.5 lc -1 title 'Exact'

plot 'out_sod_sparse.txt' u 1:(s_f($4,$2)) w p pt 6 lc "blue" title 'Sparse Large', \
     'out_sod_exact.txt' u 1:(s_f($4,$2)) w l lw 0.5 lc -1 title 'Exact'

plot 'out_sod_sparse.txt' u 1:(s_f($7,$5)) w p pt 6 lc "green" title 'Sparse Middle', \
     'out_sod_exact.txt' u 1:(s_f($4,$2)) w l lw 0.5 lc -1 title 'Exact'

plot 'out_sod_sparse.txt' u 1:(s_f($10,$8)) w p pt 6 lc "purple" title 'Sparse Low', \
     'out_sod_exact.txt' u 1:(s_f($4,$2)) w l lw 0.5 lc -1 title 'Exact'

unset multiplot
set output
exit