set term png large
set output "task3_row_runtimes.png"
set xlabel 'Number of Processors'
set ylabel 'Runtime (seconds)'
plot 'task3_row_n61.data'  with lines title "N = 61", \
    'task3_row_n121.data'  with lines title "N = 121", \
    'task3_row_n241.data'  with lines title "N = 241", \
    'task3_row_n481.data'  with lines title "N = 481", \
    'task3_row_n961.data'  with lines title "N = 961", \
    'task3_row_n1921.data'  with lines title "N = 1921"
    
set output "task3_block_runtimes.png"
plot 'task3_block_n61.data'  with lines title "N = 61", \
    'task3_block_n121.data'  with lines title "N = 121", \
    'task3_block_n241.data'  with lines title "N = 241", \
    'task3_block_n481.data'  with lines title "N = 481", \
    'task3_block_n961.data'  with lines title "N = 961", \
    'task3_block_n1921.data'  with lines title "N = 1921"    
