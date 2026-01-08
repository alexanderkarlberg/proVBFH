# gnuplot file
set term postscript enhanced color

filename='nlops-comparison.ps'
set output filename
set macros
norm=1e3
pythia6='pythia6-P11'
pythia8='pythia8-4C'
herwig='herwig-AUET2'

# filehadron=' -f tmp/lhc13-pythia8-4C-hadron.res '
# fileparton=' -f tmp/lhc13-pythia8-4C-parton.res '
# filehadrUE=' -f tmp/lhc13-pythia8-4C-hadron+UE.res '
# 
# filehadron=' -f tmp/lhc13-pythia6-P11-hadron.res '
# fileparton=' -f tmp/lhc13-pythia6-P11-parton.res '
# filehadrUE=' -f tmp/lhc13-pythia6-P11-hadron+UE.res '

# filehadron=' -f tmp/lhc13-Hgaga-pythia8-4C-hadron.res '
# fileparton=' -f tmp/lhc13-Hgaga-pythia8-4C-parton.res '
# filehadrUE=' -f tmp/lhc13-Hgaga-pythia8-4C-hadron+UE.res '

# filehadron=' -f tmp2/lhc13-Hgaga-pythia6-P11-hadron.res '
# fileparton=' -f tmp2/lhc13-Hgaga-pythia6-P11-parton.res '
# filehadrUE=' -f tmp2/lhc13-Hgaga-pythia6-P11-hadron+UE.res '

colors = '#00BFFF #50c050 #ff8000'
lc1='lc rgb word(colors,1)'
lc2='lc rgb word(colors,2)'
lc3='lc rgb word(colors,3)'

filehadr(mc)=' -f lhc13-Hgaga-powheg-NLO-'.mc.'-hadron.res '
file(mc)=' -f lhc13-Hgaga-powheg-NLO-'.mc.'-parton.res '
fileNS=' -f lhc13-Hgaga-powheg-NLO-noshower.res '
merge='< mergeidx.pl '
pasteall='< exec bash -c "paste '
me='<(mergeidx.pl '

do for [obs in 'vbfcut.HTall vbfcut.ptj1 vbfcut.ptj2 vbfcut.ptH vbfcut.yH vbfcut.yj1 vbfcut.yj2 vbfcut.dyj1j2 vbfcut.dphij1j2 vbfcut.minrap_j1j2_j2j3 vbfcut.minR_j1j2_j2j3'] {
reset
set xlabel ''
set ylabel obs.' [nb]'
set title obs
set grid
set multiplot
set origin 0.0,0.3
set size 1.0,0.7
set xtics format ""
set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.93
set bmargin at screen 0.4
plot merge.fileNS.obs            u 2:4 w lines t 'NLO',\
     merge.file(pythia6).obs     u 2:4 w lines lc rgb word(colors,1) dt 1 t 'Parton, Pythia 6 (P11)',\
     merge.file(pythia8).obs     u 2:4 w lines lc rgb word(colors,2) dt 1 t 'Parton, Pythia 8 (4C)',\
     merge.file(herwig).obs      u 2:4 w lines lc rgb word(colors,3) dt 1 t 'Parton, Herwig (AUET2)',\
     merge.filehadr(pythia6).obs u 2:4 w lines lc rgb word(colors,1) dt 2 t 'Hadron, Pythia 6 w/ P11',\
     merge.filehadr(pythia8).obs u 2:4 w lines lc rgb word(colors,2) dt 2 t 'Hadron, Pythia 8 w/ 4C',\
     merge.filehadr(herwig).obs  u 2:4 w lines lc rgb word(colors,3) dt 2 t 'Hadron, Herwig w/ AUET2'

set tmargin at screen 0.4
set bmargin at screen 0.1
set title ''
set ylabel 'ratio NLO+PS / NLO'
set xlabel 'p_T [GeV]'
set format x
set ytics 0.85,0.05,1.15
set yrange [0.80:1.20]
plot 1 w lines not,\
     pasteall.me.file(pythia6).obs.') '.me.fileNS.obs.')"'     u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 not,\
     pasteall.me.file(pythia8).obs.') '.me.fileNS.obs.')"'     u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 not,\
     pasteall.me.file(herwig).obs.') '.me.fileNS.obs.')"'      u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 not,\
     pasteall.me.filehadr(pythia6).obs.') '.me.fileNS.obs.')"' u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     pasteall.me.filehadr(pythia8).obs.') '.me.fileNS.obs.')"' u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     pasteall.me.filehadr(herwig).obs.') '.me.fileNS.obs.')"'  u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

unset multiplot
}
set output
