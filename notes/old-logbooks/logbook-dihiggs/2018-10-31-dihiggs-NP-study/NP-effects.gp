# gnuplot file
set term pdfcairo enhanced color size 16cm,8.5cm lw 2


filename='np-effects.pdf'
set output filename
set macros
norm=1e3

colors = '#00BFFF #50c050 #ff8000'
lc1='lc rgb word(colors,1)'
lc2='lc rgb word(colors,2)'
lc3='lc rgb word(colors,3)'

# hadronisation correction: hadron (no MPI) / parton
ratio_hadr=' -f lhc14-HHjj-pythia8-hadron.res -f lhc14-HHjj-pythia8-parton.res '

# Underlying event correction: hadron (MPI) / parton
ratio_mpi=' -f lhc14-HHjj-pythia8.res -f lhc14-HHjj-pythia8-hadron.res '

ratio_np=' -f lhc14-HHjj-pythia8.res -f lhc14-HHjj-pythia8-parton.res '
merge='< mergeidx.pl '
pasteall='< exec bash -c "paste '
me='<(mergeidx.pl '

rbn=' | rebin.pl -r -c 3 '

#observables: ptj1 ptj2 ptH yH dyj1j2 Mjj

# ptj1
reset
obs='vbfcut.ptj1'
print merge.ratio_hadr.obs
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
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

# ptj1
reset
obs='vbfcut.ptj2'
print merge.ratio_hadr.obs
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
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

# ptH
reset
obs='vbfcut.ptHH'
print merge.ratio_hadr.obs
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
set xlabel 'p_{t,HH} [GeV]'
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

# yH
reset
obs='vbfcut.yHH'
print merge.ratio_hadr.obs
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
set xlabel 'y_{HH} [GeV]'
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

# dyj1j2
reset
obs='vbfcut.dyj1j2'
print merge.ratio_hadr.obs
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
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

# Mjj
reset
obs='vbfcut.Mjj'
print merge.ratio_hadr.obs
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
plot merge.ratio_np.obs.rbn   u 2:($4/$9) w l lc 1 t 'UE+Hadronisation',\
     merge.ratio_hadr.obs.rbn u 2:($4/$9) w l lc 2 t 'Hadronisation',\
     merge.ratio_mpi.obs.rbn  u 2:($4/$9) w l lc 3 t 'Underlying event'

#     sqrt(-1) w lines lc rgb '#00FFFFFF' dt 1 t ' ',\
set output
