set term pdfcairo enhanced lw 2 font ',12' size 6,12 ps 0.5
set out "Num_test.pdf"

set multiplot layout 5,1

set xlabel "x"
set ylabel "rho"
set yrange [-0.1:1.1]
plot 'out_sod.txt' u 1:2 w p pt 6 lc "blue", 'out_sod_exact.txt' u 1:2 w l lc -1


set ylabel "v"

set yrange [-0.1:1.1]
plot 'out_sod.txt' u 1:3 w p pt 6 lc "blue", 'out_sod_exact.txt' u 1:3 w l lc -1

set ylabel "p"

plot 'out_sod.txt' u 1:4 w p pt 6 lc "blue", 'out_sod_exact.txt' u 1:4 w l lc -1

set ylabel "e"

e_f(r,p,u) = 1.4*r*p+0.5*r*u**2
set yrange [0:1.5]
plot 'out_sod.txt' u 1:(e_f($2,$4,$3)) w p pt 6 lc "blue", 'out_sod_exact.txt' u 1:(e_f($2,$4,$3)) w l lc -1

set ylabel "s"
set yrange [0.5:2.5]
s_f(p,r) = p/(r**1.4)
plot 'out_sod.txt' u 1:(s_f($4,$2)) w p pt 6 lc "blue", 'out_sod_exact.txt' u 1:(s_f($4,$2)) w l lc -1



unset multiplot
set out
exit