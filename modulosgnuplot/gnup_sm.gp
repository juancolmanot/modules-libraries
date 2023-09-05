terminal_p = ARG1
size_p = ARG2
datafile_p = ARG3
font_p = ARG4
output_p = ARG5
cols_p = ARG6
labels_p = ARG7
legends_p = ARG8

set_size(size_p) = sprintf("array size_p[2] = %s", size_p)
set_cols(cols_p) = sprintf("array cols_p[4] = %s", cols_p)
set_labels(labels_p) = sprintf("array labels_p[3] = %s", labels_p)
set_legends(legends_p) = sprintf("array legends_p[2] = %s", legends_p)

eval(set_size(size_p))
eval(set_cols(cols_p))
eval(set_labels(labels_p))
eval(set_legends(legends_p))

set terminal terminal_p size size_p[1], size_p[2] font font_p

set output output_p

set xlabel labels_p[1]
set ylabel labels_p[2]
set title labels_p[3]

load 'gnup_linestyles.gp'

p datafile_p u cols_p[1]:cols_p[2] t legends_p[1] w l linestyle 1, \
  datafile_p u cols_p[3]:cols_p[4] t legends_p[2] w l linestyle 2