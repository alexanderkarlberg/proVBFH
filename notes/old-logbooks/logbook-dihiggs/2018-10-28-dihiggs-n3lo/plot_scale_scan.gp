set term pdfcairo enhanced color transparent font "Times, 26" size 17cm,15cm
# set term wxt persist

set macros

LOline   = "lines lc rgb '#FFCC11' lw 3.0"
NLOline  = "lines lc rgb '#0099CC' lw 3.0"
NNLOline = "lines lc rgb '#37BC61' lw 3.0"
N3LOline = "lines lc rgb '#FF4500' lw 3.0"

pbtofb=1000
fn = 'scale_scan_27tev.dat'
fnmuf = 'scale_scan_muf_27tev.dat'
fnmur = 'scale_scan_mur_27tev.dat'
set output 'n3lo_scale_scan.pdf'
# plot cross section
reset
#set grid
set ylabel '{/Symbol s} [fb]' offset 1,0
set xlabel '{/Symbol x}' offset 0,0.3
#set label 1000 '{/*4 PRELIMINARY' at graph 0.50,0.40 center rotate by 25 tc rgb '#d0d0d0'
set label 2000 'PDF4LHC14\_nnlo\_mc' at graph 0.965,0.87 right
set label 2002 'HE-LHC 27 TeV' at graph 0.965,0.94 right
set key maxrows 2
set key bottom left Left reverse
set key width 1
set xrange [0.25:8]
set log x
set ytics 0.0,0.1,10.0
set mytics 2

set xtics add (1/8.,1/4.,1/2.,2,4,8)
set yrange [8.1:8.7]
set label 2001 '{/Symbol m}_{R,i} = {/Symbol m}_{F,i} = {/Symbol x} Q_i' at graph 0.965,0.80 right
plot fn u 1:($2*pbtofb) w @LOline   smooth mcsplines t 'LO',\
     fn u 1:($3*pbtofb) w @NLOline  smooth mcsplines t 'NLO',\
     fn u 1:($4*pbtofb) w @NNLOline smooth mcsplines t 'NNLO',\
     fn u 1:($5*pbtofb) w @N3LOline smooth mcsplines t 'N^{3}LO'

# set yrange [3.85:4.05]
# set label 2001 '{/Symbol m}_R = {/Symbol x} Q,  {/Symbol m}_F = Q' at graph 0.04,0.70 left
# plot fnmur u 1:2 w @LOline t 'LO',\
#      fnmur u 1:3 w @NLOline t 'NLO',\
#      fnmur u 1:4 w @NNLOline t 'NNLO',\
#      fnmur u 1:5 w @N3LOline t 'N3LO'

# set yrange [3.85:4.05]
# set label 2001 '{/Symbol m}_R = Q,  {/Symbol m}_F = {/Symbol x} Q' at graph 0.04,0.70 left
# plot fnmuf u 1:2 w @LOline t 'LO',\
#      fnmuf u 1:3 w @NLOline t 'NLO',\
#      fnmuf u 1:4 w @NNLOline t 'NNLO',\
#      fnmuf u 1:5 w @N3LOline t 'N3LO'

set output
