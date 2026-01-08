set output "hist_noerror.ps"
set terminal postscript color
set style data histep
# set style data lines
set key off

 #           1 : dS/dpT_j (fb/GeV)
set title "dS/dpT_j (fb/GeV)"
plot "LO/hist.1.dat" 

 #           2 : dS/dpTmax_j (fb/GeV)
set title "dS/dpTmax_j (fb/GeV)"
plot "LO/hist.2.dat" 

 #           3 : dS/dpTmin_j (fb/GeV)
set title "dS/dpTmin_j (fb/GeV)"
plot "LO/hist.3.dat" 

 #           4 : dS/dy_j (fb)
set title "dS/dy_j (fb)"
plot "LO/hist.4.dat" 

 #           5 : dS/dy_j1 (fb)
set title "dS/dy_j1 (fb)"
plot "LO/hist.5.dat" 

 #           6 : dS/dy_j2 (fb)
set title "dS/dy_j2 (fb)"
plot "LO/hist.6.dat" 

 #           7 : dS/dpTmax_l (fb/GeV)
set title "dS/dpTmax_l (fb/GeV)"
plot "LO/hist.7.dat" 

 #           8 : dS/dpTmin_l (fb/GeV)
set title "dS/dpTmin_l (fb/GeV)"
plot "LO/hist.8.dat" 

 #           9 : dS/d|eta|max_l (fb)
set title "dS/d|eta|max_l (fb)"
plot "LO/hist.9.dat" 

 #          10 : dS/d|eta|min_l (fb)
set title "dS/d|eta|min_l (fb)"
plot "LO/hist.10.dat" 

 #          11 : dS/dPhi_jj (fb)
set title "dS/dPhi_jj (fb)"
plot "LO/hist.11.dat" 

 #          12 : dS/dmHH (fb)
set title "dS/dmHH (fb)"
plot "LO/hist.12.dat" 

 #          13 : dS/dmH1 (fb)
set title "dS/dmH1 (fb)"
plot "LO/hist.13.dat" 

 #          14 : dS/dmH2 (fb)
set title "dS/dmH2 (fb)"
plot "LO/hist.14.dat" 

 #          15 : dS/dptHH (fb)
set title "dS/dptHH (fb)"
plot "LO/hist.15.dat" 

 #          16 : dS/dpTH1 (fb)
set title "dS/dpTH1 (fb)"
plot "LO/hist.16.dat" 

 #          17 : dS/dpTH2 (fb)
set title "dS/dpTH2 (fb)"
plot "LO/hist.17.dat" 

 #          18 : dS/dpTH_h (fb)
set title "dS/dpTH_h (fb)"
plot "LO/hist.18.dat" 

 #          19 : dS/dpTH_s (fb)
set title "dS/dpTH_s (fb)"
plot "LO/hist.19.dat" 

 #          20 : dS/dyHH (fb)
set title "dS/dyHH (fb)"
plot "LO/hist.20.dat" 

 #          21 : dS/dyH1 (fb)
set title "dS/dyH1 (fb)"
plot "LO/hist.21.dat" 

 #          22 : dS/dyH2 (fb)
set title "dS/dyH2 (fb)"
plot "LO/hist.22.dat" 

 #          23 : dS/dyH_h (fb)
set title "dS/dyH_h (fb)"
plot "LO/hist.23.dat" 

 #          24 : dS/dyH_s (fb)
set title "dS/dyH_s (fb)"
plot "LO/hist.24.dat" 

 #          25 : dS/dPhiHH (fb)
set title "dS/dPhiHH (fb)"
plot "LO/hist.25.dat" 

 #           1 : d2S/dy_jj dm_jj (fb/GeV)
set view map
set size 3.9/5.0,3.4/3.0
unset surface
set style data pm3d
set style function pm3d
set pm3d implicit at b
set pm3d map corners2color c1
set title "d2S/dy_jj dm_jj (fb/GeV)"
set xlabel "eta_jj"
set ylabel "m_jj"
set palette rgbformulae 30,31,32 negative
splot "LO/hist2.1.dat" 
