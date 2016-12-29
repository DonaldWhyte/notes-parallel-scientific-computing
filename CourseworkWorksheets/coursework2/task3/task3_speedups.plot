set term png large
set output "task3_row_speedups.png"
set xlabel 'Number of Processors'
set ylabel 'Speedup factor, S(p)'
plot 'task3_row_speedup_n61.data'  with lines title "N = 61", \
    'task3_row_speedup_n121.data'  with lines title "N = 121", \
    'task3_row_speedup_n241.data'  with lines title "N = 241", \
    'task3_row_speedup_n481.data'  with lines title "N = 481", \
    'task3_row_speedup_n961.data'  with lines title "N = 961", \
    'task3_row_speedup_n1921.data'  with lines title "N = 1921"
    
set output "task3_block_speedups.png"    
plot 'task3_block_speedup_n61.data'  with lines title "N = 61", \
    'task3_block_speedup_n121.data'  with lines title "N = 121", \
    'task3_block_speedup_n241.data'  with lines title "N = 241", \
    'task3_block_speedup_n481.data'  with lines title "N = 481", \
    'task3_block_speedup_n961.data'  with lines title "N = 961", \
    'task3_block_speedup_n1921.data'  with lines title "N = 1921"
