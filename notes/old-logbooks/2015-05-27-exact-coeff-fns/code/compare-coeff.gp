
set term postscript enhanced color
set output 'compare-coeff.ps'

list='F2NSplus-A F2NSplus-B F2NSminus-A F2NSminus-B FLNSplus-A FLNSminus-A F3NSplus-A F3NSplus-B F3NSminus-A F3NSminus-B F3NSWrongSignDiffplus-A'

merge='<../../../aux/mergeidx.pl -f compare-coeff.res '

set grid
set log x
set xlabel 'x'

set lmargin at screen 0.1
piece='F2NSplus-A'
do for [piece in list] {

  print piece
  set title piece
  set ylabel '(1-x) * coeff-fn'
  plot merge.piece u 1:((1-$1)*$2) t 'exact',\
       ''          u 1:((1-$1)*$3) t 'approx'
       
  set ylabel 'approx/exact'
  plot merge.piece u 1:($3/$2) t 'approx/exact',\
       1 w l lt 1 lc 0 t ''
  
}


set output
     