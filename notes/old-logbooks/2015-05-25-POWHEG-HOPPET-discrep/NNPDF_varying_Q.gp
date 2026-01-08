set term postscript color enhanced
set output 'NNPDF30_varying_Q.ps'

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 0 u (($1+$2)/2.):3 t 'bbar pos',\
     'NNPDF_from_hoppet_v_Q.top' i 0 u (($1+$2)/2.):(-$3) t 'bbar neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 1 u (($1+$2)/2.):3 t 'cbar pos',\
     'NNPDF_from_hoppet_v_Q.top' i 1 u (($1+$2)/2.):(-$3) t 'cbar neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 2 u (($1+$2)/2.):3 t 'sbar pos',\
     'NNPDF_from_hoppet_v_Q.top' i 2 u (($1+$2)/2.):(-$3) t 'sbar neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 3 u (($1+$2)/2.):3 t 'ubar pos',\
     'NNPDF_from_hoppet_v_Q.top' i 3 u (($1+$2)/2.):(-$3) t 'ubar neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 4 u (($1+$2)/2.):3 t 'dbar pos',\
     'NNPDF_from_hoppet_v_Q.top' i 4 u (($1+$2)/2.):(-$3) t 'dbar neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 5 u (($1+$2)/2.):3 t 'gluon pos',\
     'NNPDF_from_hoppet_v_Q.top' i 5 u (($1+$2)/2.):(-$3) t 'gluon neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 6 u (($1+$2)/2.):3 t 'd pos',\
     'NNPDF_from_hoppet_v_Q.top' i 6 u (($1+$2)/2.):(-$3) t 'd neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 7 u (($1+$2)/2.):3 t 'u pos',\
     'NNPDF_from_hoppet_v_Q.top' i 7 u (($1+$2)/2.):(-$3) t 'u neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 8 u (($1+$2)/2.):3 t 's pos',\
     'NNPDF_from_hoppet_v_Q.top' i 8 u (($1+$2)/2.):(-$3) t 's neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 9 u (($1+$2)/2.):3 t 'c pos',\
     'NNPDF_from_hoppet_v_Q.top' i 9 u (($1+$2)/2.):(-$3) t 'c neg'
reset     

set datafile fortran
set grid
set title 'NNPDF30\_nnlo\_as\_0118, x = 0.11 + 0.00015 * Q [0.11:0.9875]'
set log y
plot 'NNPDF_from_hoppet_v_Q.top' i 10 u (($1+$2)/2.):3 t 'b pos',\
     'NNPDF_from_hoppet_v_Q.top' i 10 u (($1+$2)/2.):(-$3) t 'b neg'
reset     

set output