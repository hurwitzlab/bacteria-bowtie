#!/usr/bin/env bash

#to be used to data munge those gff's
#put in the genome directory where your gff's are

#get rid of empty lines
egrep -v '^$' all-PATRIC.gff > allPATRIC.gff.temp
#get rid of comments
egrep -v '^#' allPATRIC.gff.temp > temp2
#only get CDS
egrep "\tCDS\t" temp2 > all-PATRIC-CDS.gff
#get the rRNA for bowtie2 (or other aligners) that support filtering by rRNA seq
egrep "\trRNA\t" temp2 > all-PATRIC-rRNA.gff
