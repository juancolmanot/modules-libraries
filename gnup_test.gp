terminal_p = ARG1
size_p = ARG2
datafile_p = ARG3
font_p = ARG4
output_p = ARG5
cols_p = ARG6
labels_p = ARG7
legends_p = ARG8
n_plots = ARG9

#print ARGV

set_size(size_p) = sprintf("array size_p[2] = %s", size_p)
set_cols(cols_p, n_plots) = sprintf("array cols_p[%d] = %s", n_plots * 2, cols_p)
set_labels(labels_p) = sprintf("array labels_p[3] = %s", labels_p)
set_legends(legends_p, n_plots) = sprintf("array legends_p[%d] = %s", n_plots * 1, legends_p)

eval(set_size(size_p))
eval(set_cols(cols_p, n_plots))
eval(set_labels(labels_p))
eval(set_legends(legends_p, n_plots))

set terminal terminal_p size size_p[1], size_p[2] font font_p

set output output_p

set xlabel labels_p[1]
set ylabel labels_p[2]
set title labels_p[3]

load '/home/juan/cursos/modulosgnuplot/gnup_linestyles.gp'

p for[i=1:n_plots] datafile_p u cols_p[i * 2 - 1]:cols_p[i * 2] t legends_p[i] linestyle i