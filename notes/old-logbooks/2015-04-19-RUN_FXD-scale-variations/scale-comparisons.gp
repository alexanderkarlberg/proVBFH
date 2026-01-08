set term postscript enhanced color
set datafile fortran
# set term wxt persist
set output 'scale-comparisons.eps'

# output from altern_method = 2
fn_max_mh = 'rebinned-max-MH.top'
fn_min_mh = 'rebinned-min-MH.top'
fn_central_mh = 'rebinned-central-MH.top'
fn_max_sqrtq1q2 = 'rebinned-max-SQRTQ1Q2.top'
fn_min_sqrtq1q2 = 'rebinned-min-SQRTQ1Q2.top'
fn_central_sqrtq1q2 = 'rebinned-central-SQRTQ1Q2.top'
fn_max_q1_q2 = 'rebinned-max-Q1_Q2.top'
fn_min_q1_q2 = 'rebinned-min-Q1_Q2.top'

#ratios
fn_ratio_max_mh = 'rebinned-ratio-max-MH.top'
fn_ratio_min_mh = 'rebinned-ratio-min-MH.top'
fn_ratio_central_mh = 'rebinned-ratio-central-MH.top'
fn_ratio_max_sqrtq1q2 = 'rebinned-ratio-max-SQRTQ1Q2.top'
fn_ratio_min_sqrtq1q2 = 'rebinned-ratio-min-SQRTQ1Q2.top'
fn_ratio_central_sqrtq1q2 = 'rebinned-ratio-central-SQRTQ1Q2.top'
fn_ratio_max_q1_q2 = 'rebinned-ratio-max-Q1_Q2.top'
fn_ratio_min_q1_q2 = 'rebinned-ratio-min-Q1_Q2.top'


# plot ptj1
reset
set grid 
set key left
set xtics format ""
set xrange [0:600]
set yrange [0.99:1.02]
set title 'ptj1 at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 7
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'ptj1 [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot ptj2
reset
set grid 
set key left
set xtics format ""
set xrange [0:600]
set yrange [0.98:1.04]
set title 'ptj2 at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 8
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'ptj2 [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''


# plot ptH
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.04]
set title 'ptH at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 9
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'ptH [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot yj1
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'yj1 at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 10
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'yj1'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''


# plot yj2
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'yj2 at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 11
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'yj2'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''


# plot yH
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'yH at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 12
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'yH'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''


# plot HT_jets
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'HT_jets at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 13
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'HT_jets [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot HT_all
reset
set grid 
set key left
set xtics format ""
set xrange [0:1000]
set yrange [0.98:1.02]
set title 'HT_all at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 14
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'HT_all [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot MHjj
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'MHjj at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 15
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'MHjj [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot Mjj
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.98:1.02]
set title 'Mjj at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 16
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'Mjj [GeV]'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot Rap(j1,j2)
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.99:1.02]
set title 'Rap(j1,j2) at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 17
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'Rap(j1,j2)'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot phi(j1,j2)
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.99:1.02]
set title 'Rap(j1,j2) at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 18
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'phi(j1,j2)'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

# plot R(j1,j2)
reset
set grid 
set key left
set xtics format ""
set xrange [*:*]
set yrange [0.99:1.02]
set title 'R(j1,j2) at NNLO, 13 TeV'
#set multiplot
set origin 0.0,0.0
set size 1.0,1.0
iindex = 19
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.15
set format x
set xlabel 'R(j1,j2)'
plot fn_ratio_max_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'max-q1_q2',\
     fn_ratio_min_q1_q2 i iindex u (($1+$2)/2):3 w line lc rgb 'blue' linewidth 2.0 ti 'min-q1_q2',\
     fn_ratio_max_sqrtq1q2 i iindex u (($1+$2)/2):3 w line lc rgb 'red' linewidth 2.0 ti 'max-sqrtq1q2',\
     fn_ratio_min_sqrtq1q2 i iindex u (($1+$2)/2):3:4 w line lc rgb 'red' linewidth 2.0 ti 'min-sqrtq1q2',\
     fn_ratio_central_sqrtq1q2 i iindex u (($1+$2)/2):3 lc rgb 'red' linewidth 2.0 ti 'cent-sqrtq1q2',\
     fn_ratio_max_mh i iindex u (($1+$2)/2):3 w line lc rgb 'magenta' linewidth 2.0 ti 'max-mh',\
     fn_ratio_min_mh i iindex u (($1+$2)/2):3:4 w line lc rgb 'magenta' linewidth 2.0 ti 'min-mh',\
     fn_ratio_central_mh i iindex u (($1+$2)/2):3 lc rgb 'magenta' linewidth 2.0 ti 'cent-mh',\
     1 lw 4 lc rgb 'blue' ti ''

set output

