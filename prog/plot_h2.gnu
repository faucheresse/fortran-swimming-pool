reset
set term pdfcairo size 64cm,36cm font 'Libertinus Serif,40'
set output "wavefunction_H2.pdf"

set grid
set xlabel "\Omega (a.u.)"
set ylabel "Population"
# set yrange [-1.4:1.4]
plot "data_H2.dat" u 1:2 w l title "1st component", "data_H2.dat" u 1:3 w l title "2nd component"
