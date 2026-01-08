set term pdfcairo enhanced transparent color font "Palatino, 14" size 15cm,15cm
set datafile fortran
set output 'unitarity-talk.pdf'

fn_TT = 'ratio-boxoff.top'
fn_BB = 'ratio-trioff.top'
fn_TB = 'ratio-interference.top'
fn_full = 'ratio-full.top'

ColTT='#e41a1c'
ColBB='#377eb8'
ColTB='#4daf4a'
ColFull='#984ea3'

#e41a1c
#377eb8
#4daf4a
#984ea3
#ff7f00
#ffff33
#a65628
#f781bf

LabTT='{/Symbol s}_{TT}'
LabBB='{/Symbol s}_{BB}'
LabTB='{/Symbol s}_{TB}'
LabFull='{/Symbol s}_{TT}+{/Symbol s}_{BB}+{/Symbol s}_{TB}'

pretitle='{/*1.2 '
tinySpace='{/*0.5 &.}'
tsize = 0.05
# plot ptj1
reset

set size 1,0.97
set origin 0,0.0

set grid front
set label 4 'LHC 13 TeV' at graph 0.05,0.756
set label 9 left "{proVBFHH v1.1.0}" at screen 0.10, screen 1.0-tsize+0.02 tc rgb "black"
set nologscale y
set title ''
set format x
set format y
set log x
#set ytics 0.50,0.2,1.5
set mytics 5
#unset label
set key left top width 2.6
unset ylabel
set xlabel font "Times, 18" 'p_{t,j_1} [GeV]' offset 35,0.5
set ylabel 'd{/Symbol s}_{nf}/d{/Symbol s}_{LO}' offset 1.
set xrange [*:*]
set yrange [*:*]
iindex=2

plot fn_TT i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColTT ti LabTT,\
     fn_BB i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColBB ti LabBB,\
     fn_TB i iindex u   (($1+$2)/2.):3 w l dt (13,5) lw 2 lc rgb ColTB ti LabTB,\
     fn_full i iindex u (($1+$2)/2.):3 w l lw 2 lc rgb ColFull ti LabFull

# plot ptj1
reset

set size 1,0.97
set origin 0,0.0

set grid front
set label 4 'LHC 13 TeV' at graph 0.05,0.75
set label 9 left "{proVBFHH v1.1.0}" at screen 0.10, screen 1.0-tsize+0.02 tc rgb "black"
set nologscale y
set title ''
set format x
set format y
set log x
set log y
#set ytics 0.50,0.2,1.5
set mytics 5
#unset label
set key left top width 2.6
unset ylabel
set xlabel font "Times, 18" 'p_{t,j_2} [GeV]' offset 35,0.5
set ylabel 'd{/Symbol s}_{nf}/d{/Symbol s}_{LO}' offset 1.
set xrange [*:*]
set yrange [*:100000]
iindex=4
plot fn_TT i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColTT ti LabTT,\
     fn_BB i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColBB ti LabBB,\
     fn_TB i iindex u   (($1+$2)/2.):3 w l dt (13,5) lw 2 lc rgb ColTB ti LabTB,\
     fn_full i iindex u (($1+$2)/2.):3 w l lw 2 lc rgb ColFull ti LabFull

# plot ptj1
reset

set size 1,0.97
set origin 0,0.0

set grid front
set label 4 'LHC 13 TeV' at graph 0.05,0.75
set label 9 left "{proVBFHH v1.1.0}" at screen 0.10, screen 1.0-tsize+0.02 tc rgb "black"
set nologscale y
set title ''
set format x
set format y
set log x
#set ytics 0.50,0.2,1.5
set mytics 5
#unset label
set key left top width 2.6
unset ylabel
set xlabel font "Times, 18" 't [GeV]' offset 35,0.5
set ylabel 'd{/Symbol s}_{nf}/d{/Symbol s}_{LO}' offset 1.
set xrange [*:4600]
set yrange [*:*]
iindex=6
plot fn_TT i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColTT ti LabTT,\
     fn_BB i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColBB ti LabBB,\
     fn_TB i iindex u   (($1+$2)/2.):3 w l dt (13,5) lw 2 lc rgb ColTB ti LabTB,\
     fn_full i iindex u (($1+$2)/2.):3 w l lw 2 lc rgb ColFull ti LabFull


# plot ptj1
reset

set size 1,0.97
set origin 0,0.0

set grid front
set label 4 'LHC 13 TeV' at graph 0.05,0.75
set label 9 left "{proVBFHH v1.1.0}" at screen 0.10, screen 1.0-tsize+0.02 tc rgb "black"
set nologscale y
set title ''
set format x
set format y
set log x
#set ytics 0.50,0.2,1.5
set mytics 5
#unset label
set key left top width 2.6
unset ylabel
set xlabel font "Times, 18" 'u [GeV]' offset 35,0.5
set ylabel 'd{/Symbol s}_{nf}/d{/Symbol s}_{LO}' offset 1.
set xrange [*:4600]
set yrange [*:*]
iindex=8
plot fn_TT i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColTT ti LabTT,\
     fn_BB i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColBB ti LabBB,\
     fn_TB i iindex u   (($1+$2)/2.):3 w l dt (13,5) lw 2 lc rgb ColTB ti LabTB,\
     fn_full i iindex u (($1+$2)/2.):3 w l lw 2 lc rgb ColFull ti LabFull


# plot ptj1
reset

set size 1,0.97
set origin 0,0.0

set grid front
set label 4 'LHC 13 TeV' at graph 0.05,0.75
set label 9 left "{proVBFHH v1.1.0}" at screen 0.10, screen 1.0-tsize+0.02 tc rgb "black"
set nologscale y
set title ''
set format x
set format y
set log x
#set ytics 0.50,0.2,1.5
set mytics 5
#unset label
set key left top width 2.6
unset ylabel
set xlabel font "Times, 18" 's [GeV]' offset 35,0.5
set ylabel 'd{/Symbol s}_{nf}/d{/Symbol s}_{LO}' offset 1.
set xrange [*:4600]
set yrange [*:*]
iindex=10

plot fn_TT i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColTT ti LabTT,\
     fn_BB i iindex u   (($1+$2)/2.):3 w l lw 2 lc rgb ColBB ti LabBB,\
     fn_TB i iindex u   (($1+$2)/2.):3 w l dt (13,5) lw 2 lc rgb ColTB ti LabTB,\
     fn_full i iindex u (($1+$2)/2.):3 w l lw 2 lc rgb ColFull ti LabFull


set output
