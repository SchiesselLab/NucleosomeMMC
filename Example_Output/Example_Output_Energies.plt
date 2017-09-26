set grid back lw 6 lt 0 lc rgb "#AAAAAA"
set border lw 6
set bmargin 2.4
set tmargin 2

set encoding utf8

set xlabel "Sweep number"  font "Helvetica-Bold, 20" offset 0,1.2
set ylabel "System Energy (kT)"  font "Helvetica-Bold, 20" offset 2,0

set xtics font "Helvetica-Bold, 18" offset 0,0.5
set ytics font "Helvetica-Bold, 18" offset 0.8,0

set xrange [0:100]

set key above

set term pdf enhanced font "Helvetica-Bold, 20" dashed size 6,6
set output 'Sample_Output_Energies.pdf'

plot 'Example_Output_Energies.txt' u 1:2 w l lw 6 lt 2 lc rgb 'forest-green' notitle
     
unset output

