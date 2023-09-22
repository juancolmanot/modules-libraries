terminal_p = ARG1
size_p = ARG2
datafile_p = ARG3
font_p = ARG4
output_p = ARG5
cols_p = ARG6
labels_p = ARG7
legend_p = ARG8

set_size(size_p) = sprintf("array size_p[2] = %s", size_p)
set_cols(cols_p) = sprintf("array cols_p[2] = %s", cols_p)
set_labels(labels_p) = sprintf("array labels_p[3] = %s", labels_p)

eval(set_size(size_p))
eval(set_cols(cols_p))
eval(set_labels(labels_p))

set terminal terminal_p size size_p[1], size_p[2] font font_p

set output output_p


set xlabel labels_p[1]
set ylabel labels_p[2]
set title labels_p[3]

load 'gnup_linestyles.gp'

p datafile_p u cols_p[1]:cols_p[2] t legend_p w l linestyle 1