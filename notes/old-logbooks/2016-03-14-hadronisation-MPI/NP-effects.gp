# gnuplot file
set term pdfcairo enhanced color size 16cm,8.5cm lw 2


filename='np-effects.pdf'
set output filename
set macros
norm=1e3
pythia6='pythia6-P11'
pythia8='pythia8-4C'
herwig='herwig-AUET2'

colors = '#00BFFF #50c050 #ff8000'
lc1='lc rgb word(colors,1)'
lc2='lc rgb word(colors,2)'
lc3='lc rgb word(colors,3)'

fileNS=' -f lhc13-Hgaga-powheg-NLO-noshower.res '

# hadronisation correction: hadron (no MPI) / parton
ratio_hadr(mc)=' -f lhc13-Hgaga-powheg-NLO-'.mc.'-hadron.res -f lhc13-Hgaga-powheg-NLO-'.mc.'-parton.res '

# Underlying event correction: hadron (MPI) / parton
ratio_mpi(mc)=' -f lhc13-Hgaga-powheg-NLO-'.mc.'-hadron-UE.res -f lhc13-Hgaga-powheg-NLO-'.mc.'-hadron.res '

ratio_np(mc)=' -f lhc13-Hgaga-powheg-NLO-'.mc.'-hadron-UE.res -f lhc13-Hgaga-powheg-NLO-'.mc.'-parton.res '
merge='< mergeidx.pl '
pasteall='< exec bash -c "paste '
me='<(mergeidx.pl '

#observables: ptj1 ptj2 ptH yH dyj1j2 Mjj

# ptj1
reset
obs='vbfcut.ptj1'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.94 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [20:200]
set yrange [0.9:1.1]
set xlabel 'p_{t,j_1} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

# ptj1
reset
obs='vbfcut.ptj2'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.94 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [20:140]
set yrange [0.9:1.1]
set xlabel 'p_{t,j_2} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

# ptH
reset
obs='vbfcut.ptH'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.94 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [0:250]
set yrange [0.9:1.1]
set xlabel 'p_{t,H} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

# yH
reset
obs='vbfcut.yH'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.065 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [0:3.5]
set yrange [0.9:1.1]
set xlabel 'y_{H} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

# dyj1j2
reset
obs='vbfcut.dyj1j2'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.94 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [4.5:9]
set yrange [0.9:1.1]
set xlabel '{/Symbol D}y_{j_1,j_2} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

# Mjj
reset
obs='vbfcut.Mjj'
print merge.ratio_hadr(pythia6).obs
set ylabel 'Non-perturbative corrections'
set label 2 'VBF cuts' at graph 0.98,0.94 right
#set title obs
set grid
set key maxrows 3
set key opaque
set key bottom left
set ytics 0.85,0.05,1.15
set xrange [600:2400]
set yrange [0.9:1.1]
set xlabel 'm_{j_1,j_2} [GeV]'
plot sqrt(-1) w lines lc rgb 'black' dt 1 t 'Hadronisation + UE',\
     sqrt(-1) w lines lc rgb 'black' dt 4 t 'Hadronisation',\
     sqrt(-1) w lines lc rgb 'black' dt 2 t 'Underlying event',\
     merge.ratio_np(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 1 t 'Pythia 6 (P11)',\
     merge.ratio_np(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 1 t 'Pythia 8 (4C)',\
     merge.ratio_np(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 1 t 'Herwig 6 (AUET2)',\
     merge.ratio_hadr(pythia6).obs u 2:($4/$9) w lines lc rgb word(colors,1) dt 4 not,\
     merge.ratio_hadr(pythia8).obs u 2:($4/$9) w lines lc rgb word(colors,2) dt 4 not,\
     merge.ratio_hadr(herwig).obs  u 2:($4/$9) w lines lc rgb word(colors,3) dt 4 not,\
     merge.ratio_mpi(pythia6).obs  u 2:($4/$9) w lines lc rgb word(colors,1) dt 2 not,\
     merge.ratio_mpi(pythia8).obs  u 2:($4/$9) w lines lc rgb word(colors,2) dt 2 not,\
     merge.ratio_mpi(herwig).obs   u 2:($4/$9) w lines lc rgb word(colors,3) dt 2 not

#     sqrt(-1) w lines lc rgb '#00FFFFFF' dt 1 t ' ',\
set output
