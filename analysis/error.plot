set title "Maximum error over domain [0,1]x[0.1,1]"
set ylabel "Log of max error"
set xlabel "J"
set key right

set xrange [9:201]

set style function linespoints
set terminal pdf font "Veranda, 14"
set output "error_plot.pdf"
plot '1.txt' with linespoints title "θ=0, µ=0.5", \
    '2.txt' with linespoints  title "θ=1, µ=5",     \
    '3.txt' with linespoints  title "θ=0.5, ν=1/20" 
