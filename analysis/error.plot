set ylabel "Log of max error"
set key right
set style function linespoints
set logscale x
set terminal pdf font "Veranda, 14"


set output "error_vs_J.pdf"
set title "Maximum error over domain [0,1] x [0.1,1] vs spatial grid dimension"
set xlabel "J"
set xtics (10,20,30,40,50,60,70,80)
set xrange [9:81]
plot '1.txt' using 1:3 with linespoints  title "θ=0, µ=0.5",   \
     '2.txt' using 1:3 with linespoints  title "θ=1, µ=5",     \
     '3.txt' using 1:3 with linespoints  title "θ=0.5, ν=1/20" 

set output "error_vs_grid.pdf"
set title "Maximum error over domain [0,1] x [0.1,1] vs number of grid points"
set xlabel "Number of grid points"
set xrange [100:1000000]
set xtics (100,1000,10000,100000,1000000)
plot '1.txt' using 2:3 with linespoints  title "θ=0, µ=0.5",   \
     '2.txt' using 2:3 with linespoints  title "θ=1, µ=5",     \
     '3.txt' using 2:3 with linespoints  title "θ=0.5, ν=1/20" 
