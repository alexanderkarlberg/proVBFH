set term pdfcairo enhanced color transparent font "Times, 26" size 16cm,12cm
# set term wxt persist

set macros

LOline   = "lines lc rgb '#FFCC11' lw 3.0"
NLOline  = "lines lc rgb '#0099CC' lw 3.0"
NNLOline = "lines lc rgb '#37BC61' lw 3.0"
N3LOline = "lines lc rgb '#FF4500' lw 3.0"

fn = 'n3lo_scale-scan_Q.dat'
fnmuf = 'n3lo_scale-scan-muf_Q.dat'
fnmur = 'n3lo_scale-scan-mur_Q.dat'
set output 'n3lo_scale-scan_Q.pdf'
# plot cross section
reset
#set grid
set title '{/Symbol s} [pb]' offset 0,-0.5
set xlabel '{/Symbol x}' offset 0,0.3
#set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'PDF4LHC14\_nnlo\_mc' at graph 0.04,0.80 left
set label 2002 'LHC 13 TeV' at graph 0.04,0.90 left
set key maxrows 2
set key bottom left Left reverse
set key width 3
set xrange [0.25:4]
set log x
set ytics 3,0.05,4.2
set mytics 5

set xtics add (1/8.,1/4.,1/2.,2,4,8)
set yrange [3.85:4.05]
set label 2001 '{/Symbol m}_R = {/Symbol m}_F = {/Symbol x} Q' at graph 0.04,0.70 left
plot fn u 1:2 w @LOline t 'LO',\
     fn u 1:3 w @NLOline t 'NLO',\
     fn u 1:4 w @NNLOline t 'NNLO',\
     fn u 1:5 w @N3LOline t 'N3LO'

set yrange [3.85:4.05]
set label 2001 '{/Symbol m}_R = {/Symbol x} Q,  {/Symbol m}_F = Q' at graph 0.04,0.70 left
plot fnmur u 1:2 w @LOline t 'LO',\
     fnmur u 1:3 w @NLOline t 'NLO',\
     fnmur u 1:4 w @NNLOline t 'NNLO',\
     fnmur u 1:5 w @N3LOline t 'N3LO'

set yrange [3.85:4.05]
set label 2001 '{/Symbol m}_R = Q,  {/Symbol m}_F = {/Symbol x} Q' at graph 0.04,0.70 left
plot fnmuf u 1:2 w @LOline t 'LO',\
     fnmuf u 1:3 w @NLOline t 'NLO',\
     fnmuf u 1:4 w @NNLOline t 'NNLO',\
     fnmuf u 1:5 w @N3LOline t 'N3LO'

set output
