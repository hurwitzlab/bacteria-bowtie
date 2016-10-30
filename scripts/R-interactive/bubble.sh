#!/usr/bin/env bash

# Wrapper for bubble.pl
# So that it autoconverts
# to a .pdf for use in illustrator
# requires inkscape obviously

if [[ $# = 0 ]]; then
  echo "First argument is your .csv file"
  echo "Second argument is the NAME of you .pdf file you want"
  echo "WITHOUT the extension"
fi

perl ./bubble.pl -o $2.svg $1
sleep 1
inkscape $2.svg --export-pdf=$2.pdf
