set logscale xy
set term png large
set output "task1_plot.png"
set xlabel 'log(n)'
set ylabel 'log(|error|)'
plot 'task1.data' with lines
