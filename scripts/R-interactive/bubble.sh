#!/usr/bin/env bash

# Wrapper for bubble.pl
# So that it autoconverts
# to a .pdf for use in illustrator
# requires inkscape obviously

if [[ $# = 0 ]]; then
    echo "First argument is your .tab file"
    echo "Second argument is the NAME of you .pdf file you want"
    echo "WITHOUT the extension"
    exit 1
fi

perl ./bubble.pl -r -s -o $2.svg $1
sleep 1
inkscape $2.svg --export-pdf=$2.pdf
sleep 1
rm $2.svg
