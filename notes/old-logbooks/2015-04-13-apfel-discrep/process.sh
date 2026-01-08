#!/bin/zsh
# script to convert Apfel's .C files into .
# It assumes a format with 100 lines...

if [[ $# -ne 1 ]]; then
    echo "Usage: ./process.sh rootfile.C"
    exit -1
fi

infile=$1
outfile=$1:s/\.C/.dat/
echo Sending output to $outfile

echo "# "$* >! $outfile
echo "# i  x  result" >> $outfile
cat $1 | \
    sed -n '/Graph_Graph1/q;p' | \
    grep 'SetPoint(' | \
    sed -e 's/.*(//' -e 's/).*//' -e 's/,/ /g' >> $outfile
