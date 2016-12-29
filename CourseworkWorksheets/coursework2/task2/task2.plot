set logscale xy
set term png large
set output "task2_plot.png"
set xlabel 'log(n)'
set ylabel 'log(|error|)'
plot 'task2.data' with lines
