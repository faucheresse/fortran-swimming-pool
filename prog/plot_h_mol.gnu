reset
set term pdfcairo size 64cm,36cm font 'Libertinus Serif,40'
set output "wavefunction_Hmol.pdf"

set grid
set xlabel "L (a.u.)"
set ylabel "Population"
# set yrange [-1.4:1.4]
plot "data_Hmol.dat" u 1:2 w l title "1st component", "data_Hmol.dat" u 1:3 w l title "2nd component", "data_Hmol.dat" u 1:4 w l title "3rd component"
