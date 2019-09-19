set key bottom center
set style function linespoints
set terminal pdf font "Veranda, 14"
set xrange [0:1]


set output "solution1.pdf"
set title "Forward Euler scheme (µ=0.5) vs exact solution at t = 0.5"
plot '1_sol.txt' using 1:2 with points title "Numerical approximation",   \
     '1_sol.txt' using 1:3 with lines  title "Exact solution" 

set output "solution2.pdf"
set title "Backward Euler scheme (µ=5) vs exact solution at t = 0.5"
plot '2_sol.txt' using 1:2 with points title "Numerical approximation",   \
     '2_sol.txt' using 1:3 with lines  title "Exact solution" 

set output "solution3.pdf"
set title "Crank-Nicolson scheme (ν=0.05) vs exact solution at t = 0.5"
plot '3_sol.txt' using 1:2 with points title "Numerical approximation",   \
     '3_sol.txt' using 1:3 with lines  title "Exact solution" 
