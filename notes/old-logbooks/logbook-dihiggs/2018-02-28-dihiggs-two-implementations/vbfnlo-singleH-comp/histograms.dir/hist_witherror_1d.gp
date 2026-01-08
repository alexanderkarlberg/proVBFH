set output "hist_witherror_1d.ps"
set terminal postscript color
set style data histep
# set style data lines
set key off

 #           1 : dS/dpT_j (fb/GeV)
set title "dS/dpT_j (fb/GeV)"
plot "LO/hist.1.dat", "LO/hist.1.dat" with yerrorbars lt 1 pt 0

 #           2 : dS/dpTmax_j (fb/GeV)
set title "dS/dpTmax_j (fb/GeV)"
plot "LO/hist.2.dat", "LO/hist.2.dat" with yerrorbars lt 1 pt 0

 #           3 : dS/dpTmin_j (fb/GeV)
set title "dS/dpTmin_j (fb/GeV)"
plot "LO/hist.3.dat", "LO/hist.3.dat" with yerrorbars lt 1 pt 0

 #           4 : dS/dy_j (fb)
set title "dS/dy_j (fb)"
plot "LO/hist.4.dat", "LO/hist.4.dat" with yerrorbars lt 1 pt 0

 #           5 : dS/dy_j1 (fb)
set title "dS/dy_j1 (fb)"
plot "LO/hist.5.dat", "LO/hist.5.dat" with yerrorbars lt 1 pt 0

 #           6 : dS/dy_j2 (fb)
set title "dS/dy_j2 (fb)"
plot "LO/hist.6.dat", "LO/hist.6.dat" with yerrorbars lt 1 pt 0

 #           7 : dS/dpTmax_l (fb/GeV)
set title "dS/dpTmax_l (fb/GeV)"
plot "LO/hist.7.dat", "LO/hist.7.dat" with yerrorbars lt 1 pt 0

 #           8 : dS/dpTmin_l (fb/GeV)
set title "dS/dpTmin_l (fb/GeV)"
plot "LO/hist.8.dat", "LO/hist.8.dat" with yerrorbars lt 1 pt 0

 #           9 : dS/d|eta|max_l (fb)
set title "dS/d|eta|max_l (fb)"
plot "LO/hist.9.dat", "LO/hist.9.dat" with yerrorbars lt 1 pt 0

 #          10 : dS/d|eta|min_l (fb)
set title "dS/d|eta|min_l (fb)"
plot "LO/hist.10.dat", "LO/hist.10.dat" with yerrorbars lt 1 pt 0

 #          11 : dS/dPhi_jj (fb)
set title "dS/dPhi_jj (fb)"
plot "LO/hist.11.dat", "LO/hist.11.dat" with yerrorbars lt 1 pt 0

 #          12 : dS/dmHH (fb)
set title "dS/dmHH (fb)"
plot "LO/hist.12.dat", "LO/hist.12.dat" with yerrorbars lt 1 pt 0

 #          13 : dS/dmH1 (fb)
set title "dS/dmH1 (fb)"
plot "LO/hist.13.dat", "LO/hist.13.dat" with yerrorbars lt 1 pt 0

 #          14 : dS/dmH2 (fb)
set title "dS/dmH2 (fb)"
plot "LO/hist.14.dat", "LO/hist.14.dat" with yerrorbars lt 1 pt 0

 #          15 : dS/dptHH (fb)
set title "dS/dptHH (fb)"
plot "LO/hist.15.dat", "LO/hist.15.dat" with yerrorbars lt 1 pt 0

 #          16 : dS/dpTH1 (fb)
set title "dS/dpTH1 (fb)"
plot "LO/hist.16.dat", "LO/hist.16.dat" with yerrorbars lt 1 pt 0

 #          17 : dS/dpTH2 (fb)
set title "dS/dpTH2 (fb)"
plot "LO/hist.17.dat", "LO/hist.17.dat" with yerrorbars lt 1 pt 0

 #          18 : dS/dpTH_h (fb)
set title "dS/dpTH_h (fb)"
plot "LO/hist.18.dat", "LO/hist.18.dat" with yerrorbars lt 1 pt 0

 #          19 : dS/dpTH_s (fb)
set title "dS/dpTH_s (fb)"
plot "LO/hist.19.dat", "LO/hist.19.dat" with yerrorbars lt 1 pt 0

 #          20 : dS/dyHH (fb)
set title "dS/dyHH (fb)"
plot "LO/hist.20.dat", "LO/hist.20.dat" with yerrorbars lt 1 pt 0

 #          21 : dS/dyH1 (fb)
set title "dS/dyH1 (fb)"
plot "LO/hist.21.dat", "LO/hist.21.dat" with yerrorbars lt 1 pt 0

 #          22 : dS/dyH2 (fb)
set title "dS/dyH2 (fb)"
plot "LO/hist.22.dat", "LO/hist.22.dat" with yerrorbars lt 1 pt 0

 #          23 : dS/dyH_h (fb)
set title "dS/dyH_h (fb)"
plot "LO/hist.23.dat", "LO/hist.23.dat" with yerrorbars lt 1 pt 0

 #          24 : dS/dyH_s (fb)
set title "dS/dyH_s (fb)"
plot "LO/hist.24.dat", "LO/hist.24.dat" with yerrorbars lt 1 pt 0

 #          25 : dS/dPhiHH (fb)
set title "dS/dPhiHH (fb)"
plot "LO/hist.25.dat", "LO/hist.25.dat" with yerrorbars lt 1 pt 0
