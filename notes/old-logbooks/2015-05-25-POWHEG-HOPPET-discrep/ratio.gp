set term postscript enhanced color 
set datafile fortran
set output 'ratio.ps'

set yrange [0.95:1.05]
set ylabel 'HOPPET/POWHEG'
set xlabel 'pt_{j_2} [GeV]'
iindex = 8

plot 1 lw 2 lc rgb 'black' title '',\
     'ratio_MSTW_HOPPET_POWHEG.top' i iindex u (($1+$2)/2.):3:4 w e lw 2 lc rgb 'red' ti 'MSTW',\
     'ratio_CT10_HOPPET_POWHEG.top' i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'CT10',\
     'ratio_NNPDF_HOPPET_POWHEG.top' i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'NNPDF',\
     'ratio_NNPDF_HOPPET_POWHEG_new.top' i iindex u (($1+$2)/2.):3:4 w e lw 2 ti 'NNPDF w neg'


set output

