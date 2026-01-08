# gnuplot file


reset
set sty data linespoints

set term postscript enhanced colour solid lw 1

set macros

badfactor=0.995
LOWp='($15-$2)/$2'
NLOWp='($15+$17-$2-$4)/$2'
NNLOWp='($15+$17+$19-$2-$4-$6)/$2'
N3LOWp='($15+$17+$19+$21-$2-$4-$6-$8)/$2'
N3LObadWp=sprintf('($15+$17+$19+$21-$2-$4-$6-%f*$8)/$2',badfactor)
LOWm='($16-$3)/$3'
NLOWm='($16+$18-$3-$5)/$3'
NNLOWm='($16+$18+$20-$3-$5-$7)/$3'
N3LOWm='($16+$18+$20+$22-$3-$5-$7-$9)/$3'
N3LObadWm=sprintf('($16+$18+$20+$22-$3-$5-$7-%f*$9)/$3',badfactor)
LOZ='($23-$10)/$10'
NLOZ='($23+$24-$10-$11)/$10'
NNLOZ='($23+$24+$25-$10-$11-$12)/$10'
N3LOZ='($23+$24+$25+$26-$10-$11-$12-$13)/$10'
N3LObadZ=sprintf('($23+$24+$25+$26-$10-$11-$12-%f*$13)/$10',badfactor)

do for [scaletyp in 'scQ scMH scMIX'] {
do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Wp_'.scaletyp.'.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='hoppet-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_'.scaletyp.'.dat'
filepair(as)='<paste str_fct_as'.as.'/'.base(xR1,xF1).' str_fct_as'.as.'/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'
print filepair('0.1024')
# pow=1
# set title 'Check at LO'
# set label 1 LOWp
# plot filepair('0.1024') u 1:(  1**pow*@LOWp),\
#      filepair('0.0512') u 1:(  2**pow*@LOWp),\
#      filepair('0.0256') u 1:(  4**pow*@LOWp),\
#      filepair('0.0128') u 1:(  8**pow*@LOWp),\
#      filepair('0.0064') u 1:( 16**pow*@LOWp),\
#      filepair('0.0032') u 1:( 32**pow*@LOWp),\
#      filepair('0.0016') u 1:( 64**pow*@LOWp),\
#      filepair('0.0008') u 1:(128**pow*@LOWp)


# pow=2
# set title 'Check at NLO'
# set label 1 NLOWp
# plot filepair('0.1024') u 1:(  1**pow*@NLOWp),\
#      filepair('0.0512') u 1:(  2**pow*@NLOWp),\
#      filepair('0.0256') u 1:(  4**pow*@NLOWp),\
#      filepair('0.0128') u 1:(  8**pow*@NLOWp),\
#      filepair('0.0064') u 1:( 16**pow*@NLOWp),\
#      filepair('0.0032') u 1:( 32**pow*@NLOWp),\
#      filepair('0.0016') u 1:( 64**pow*@NLOWp),\
#      filepair('0.0008') u 1:(128**pow*@NLOWp)


# pow=3
# set title 'Check at NNLO'
# set label 1 NNLOWp
# plot filepair('0.1024') u 1:(  1**pow*@NNLOWp),\
#      filepair('0.0512') u 1:(  2**pow*@NNLOWp),\
#      filepair('0.0256') u 1:(  4**pow*@NNLOWp),\
#      filepair('0.0128') u 1:(  8**pow*@NNLOWp),\
#      filepair('0.0064') u 1:( 16**pow*@NNLOWp),\
#      filepair('0.0032') u 1:( 32**pow*@NNLOWp),\
#      filepair('0.0016') u 1:( 64**pow*@NNLOWp),\
#      filepair('0.0008') u 1:(128**pow*@NNLOWp)

pow=4
set title 'Check at N3LO'
set label 1 N3LOWp
plot filepair('0.1024') u 1:(  1**pow*@N3LOWp),\
     filepair('0.0512') u 1:(  2**pow*@N3LOWp),\
     filepair('0.0256') u 1:(  4**pow*@N3LOWp),\
     filepair('0.0128') u 1:(  8**pow*@N3LOWp),\
     filepair('0.0064') u 1:( 16**pow*@N3LOWp),\
     filepair('0.0032') u 1:( 32**pow*@N3LOWp),\
     filepair('0.0016') u 1:( 64**pow*@N3LOWp),\
     filepair('0.0008') u 1:(128**pow*@N3LOWp)

set title sprintf('Check of artificial bad N3LO (with N3LO piece scaled by %f)',badfactor)
set label 1 N3LObadWp
plot filepair('0.1024') u 1:(  1**pow*@N3LObadWp),\
     filepair('0.0512') u 1:(  2**pow*@N3LObadWp),\
     filepair('0.0256') u 1:(  4**pow*@N3LObadWp),\
     filepair('0.0128') u 1:(  8**pow*@N3LObadWp),\
     filepair('0.0064') u 1:( 16**pow*@N3LObadWp),\
     filepair('0.0032') u 1:( 32**pow*@N3LObadWp),\
     filepair('0.0016') u 1:( 64**pow*@N3LObadWp),\
     filepair('0.0008') u 1:(128**pow*@N3LObadWp)
}
set output
}

do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Wm_'.scaletyp.'.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='hoppet-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_'.scaletyp.'.dat'
filepair(as)='<paste str_fct_as'.as.'/'.base(xR1,xF1).' str_fct_as'.as.'/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'

pow=1
# set title 'Check at LO'
# set label 1 LOWm
# plot filepair('0.1024') u 1:(  1**pow*@LOWm),\
#      filepair('0.0512') u 1:(  2**pow*@LOWm),\
#      filepair('0.0256') u 1:(  4**pow*@LOWm),\
#      filepair('0.0128') u 1:(  8**pow*@LOWm),\
#      filepair('0.0064') u 1:( 16**pow*@LOWm),\
#      filepair('0.0032') u 1:( 32**pow*@LOWm),\
#      filepair('0.0016') u 1:( 64**pow*@LOWm),\
#      filepair('0.0008') u 1:(128**pow*@LOWm)


# pow=2
# set title 'Check at NLO'
# set label 1 NLOWm
# plot filepair('0.1024') u 1:(  1**pow*@NLOWm),\
#      filepair('0.0512') u 1:(  2**pow*@NLOWm),\
#      filepair('0.0256') u 1:(  4**pow*@NLOWm),\
#      filepair('0.0128') u 1:(  8**pow*@NLOWm),\
#      filepair('0.0064') u 1:( 16**pow*@NLOWm),\
#      filepair('0.0032') u 1:( 32**pow*@NLOWm),\
#      filepair('0.0016') u 1:( 64**pow*@NLOWm),\
#      filepair('0.0008') u 1:(128**pow*@NLOWm)


# pow=3
# set title 'Check at NNLO'
# set label 1 NNLOWm
# plot filepair('0.1024') u 1:(  1**pow*@NNLOWm),\
#      filepair('0.0512') u 1:(  2**pow*@NNLOWm),\
#      filepair('0.0256') u 1:(  4**pow*@NNLOWm),\
#      filepair('0.0128') u 1:(  8**pow*@NNLOWm),\
#      filepair('0.0064') u 1:( 16**pow*@NNLOWm),\
#      filepair('0.0032') u 1:( 32**pow*@NNLOWm),\
#      filepair('0.0016') u 1:( 64**pow*@NNLOWm),\
#      filepair('0.0008') u 1:(128**pow*@NNLOWm)


pow=4
set title 'Check at N3LO'
set label 1 N3LOWm
plot filepair('0.1024') u 1:(  1**pow*@N3LOWm),\
     filepair('0.0512') u 1:(  2**pow*@N3LOWm),\
     filepair('0.0256') u 1:(  4**pow*@N3LOWm),\
     filepair('0.0128') u 1:(  8**pow*@N3LOWm),\
     filepair('0.0064') u 1:( 16**pow*@N3LOWm),\
     filepair('0.0032') u 1:( 32**pow*@N3LOWm),\
     filepair('0.0016') u 1:( 64**pow*@N3LOWm),\
     filepair('0.0008') u 1:(128**pow*@N3LOWm)

set title sprintf('Check of artificial bad N3LO (with N3LO piece scaled by %f)',badfactor)
set label 1 N3LObadWm
plot filepair('0.1024') u 1:(  1**pow*@N3LObadWm),\
     filepair('0.0512') u 1:(  2**pow*@N3LObadWm),\
     filepair('0.0256') u 1:(  4**pow*@N3LObadWm),\
     filepair('0.0128') u 1:(  8**pow*@N3LObadWm),\
     filepair('0.0064') u 1:( 16**pow*@N3LObadWm),\
     filepair('0.0032') u 1:( 32**pow*@N3LObadWm),\
     filepair('0.0016') u 1:( 64**pow*@N3LObadWm),\
     filepair('0.0008') u 1:(128**pow*@N3LObadWm)
}
set output
}

do for [structf in '1 2 3']{
filename='check-scale-dep-f'.structf.'Z_'.scaletyp.'.ps'
set output filename

# loop over scale choices in format
# ab_cd => xmuR1=a, xmuF1=b, xmuR2=c, xmuF2=d
scalechoices='11_14 11_41 11_43 14_41 14_43 41_43'

do for [scale in scalechoices] {
xR1 = scale[1:1]
xF1 = scale[2:2]
xR2 = scale[4:4]
xF2 = scale[5:5]
print 'xmuR1 = '.xR1.', xmuF1 = '.xF1.', xmuR2 = '.xR2.', xmuF2 = '.xF2
base(r,f)='hoppet-F'.structf.'WpmZ_mur'.r.'.0_muf'.f.'.0_'.scaletyp.'.dat'
filepair(as)='<paste str_fct_as'.as.'/'.base(xR1,xF1).' str_fct_as'.as.'/'.base(xR2,xF2)
set log x

set xrange [0.01:0.8]
set xlabel 'x'
set ylabel '(scale-dependence/LO-result) * 1/{/Symbol a}^{p+1} [arb. units]'


set label 1 at graph 0.05,0.05
set label 2 at graph 0.05,0.10 '{/*0.7 '.system('echo "'.filepair('0.1024').'" | sed "s/_/\\\_/g"')
set label 3 at graph 0.05,0.15 'As {/Symbol a} -> 0, curves should tend to overlap'
set label 4 at graph 0.05,0.20 '(xmuR, xmuF) = ('.xR1.','.xF1.'), ('.xR2.','.xF2.')'

# pow=1
# set title 'Check at LO'
# set label 1 LOZ
# plot filepair('0.1024') u 1:(  1**pow*@LOZ),\
#      filepair('0.0512') u 1:(  2**pow*@LOZ),\
#      filepair('0.0256') u 1:(  4**pow*@LOZ),\
#      filepair('0.0128') u 1:(  8**pow*@LOZ),\
#      filepair('0.0064') u 1:( 16**pow*@LOZ),\
#      filepair('0.0032') u 1:( 32**pow*@LOZ),\
#      filepair('0.0016') u 1:( 64**pow*@LOZ),\
#      filepair('0.0008') u 1:(128**pow*@LOZ)


# pow=2
# set title 'Check at NLO'
# set label 1 NLOZ
# plot filepair('0.1024') u 1:(  1**pow*@NLOZ),\
#      filepair('0.0512') u 1:(  2**pow*@NLOZ),\
#      filepair('0.0256') u 1:(  4**pow*@NLOZ),\
#      filepair('0.0128') u 1:(  8**pow*@NLOZ),\
#      filepair('0.0064') u 1:( 16**pow*@NLOZ),\
#      filepair('0.0032') u 1:( 32**pow*@NLOZ),\
#      filepair('0.0016') u 1:( 64**pow*@NLOZ),\
#      filepair('0.0008') u 1:(128**pow*@NLOZ)


# pow=3
# set title 'Check at NNLO'
# set label 1 NNLOZ
# plot filepair('0.1024') u 1:(  1**pow*@NNLOZ),\
#      filepair('0.0512') u 1:(  2**pow*@NNLOZ),\
#      filepair('0.0256') u 1:(  4**pow*@NNLOZ),\
#      filepair('0.0128') u 1:(  8**pow*@NNLOZ),\
#      filepair('0.0064') u 1:( 16**pow*@NNLOZ),\
#      filepair('0.0032') u 1:( 32**pow*@NNLOZ),\
#      filepair('0.0016') u 1:( 64**pow*@NNLOZ),\
#      filepair('0.0008') u 1:(128**pow*@NNLOZ)

pow=4
set title 'Check at N3LO'
set label 1 N3LOZ
plot filepair('0.1024') u 1:(  1**pow*@N3LOZ),\
     filepair('0.0512') u 1:(  2**pow*@N3LOZ),\
     filepair('0.0256') u 1:(  4**pow*@N3LOZ),\
     filepair('0.0128') u 1:(  8**pow*@N3LOZ),\
     filepair('0.0064') u 1:( 16**pow*@N3LOZ),\
     filepair('0.0032') u 1:( 32**pow*@N3LOZ),\
     filepair('0.0016') u 1:( 64**pow*@N3LOZ),\
     filepair('0.0008') u 1:(128**pow*@N3LOZ)

set title sprintf('Check of artificial bad N3LO (with N3LO piece scaled by %f)',badfactor)
set label 1 N3LObadZ
plot filepair('0.1024') u 1:(  1**pow*@N3LObadZ),\
     filepair('0.0512') u 1:(  2**pow*@N3LObadZ),\
     filepair('0.0256') u 1:(  4**pow*@N3LObadZ),\
     filepair('0.0128') u 1:(  8**pow*@N3LObadZ),\
     filepair('0.0064') u 1:( 16**pow*@N3LObadZ),\
     filepair('0.0032') u 1:( 32**pow*@N3LObadZ),\
     filepair('0.0016') u 1:( 64**pow*@N3LObadZ),\
     filepair('0.0008') u 1:(128**pow*@N3LObadZ)
}
set output
}
}
