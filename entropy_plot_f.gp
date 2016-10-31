set term pdf enhanced
set output "BNEWlinear2.pdf"
set xlabel "t" font "Times-New Roman,20"
set ylabel "S(t)" font "Times-New Roman,20"
set termoption enhanced
set key bottom right
set xtics font "Times-New Roman,20"
set ytics font "Times-New Roman,20"
set key font "Times-New Roman,17"
set size ratio 0.7

plot "ISING_test3.txt" u 1:2 w l lt rgb "orange" lw 5.3 title " N = 34 ",\
"ISING_CH50spinsto1.txt" u 1:2 w l lt rgb "blue" lw 5.3 title " N = 50 "
