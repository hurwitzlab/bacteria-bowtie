#echo STARTED $(date)
#touch RNA_counts
#echo "RNA_1" >> RNA_counts
#find ./ -iname \*RNA_1\* -print0 | xargs -0 -I file egrep '^@.*' file \
#   | wc -l >> RNA_counts
#echo "Counting reads in RNA_2"
#echo "RNA_2" >> RNA_counts
#find ./ -iname \*RNA_2\* -print0 | xargs -0 -I file egrep '^@.*' file \
#    | wc -l >> RNA_counts
#echo "Counting reads in RNA_3"
#echo "RNA_3" >> RNA_counts
#find ./ -iname \*RNA_3\* -print0 | xargs -0 -I file egrep '^@.*' file \
#    | wc -l >> RNA_counts
#echo "Counting reads in RNA_4"
#echo "RNA_4" >> RNA_counts
#find ./ -iname \*RNA_4\* -print0 | xargs -0 -I file egrep '^@.*' file \
#    | wc -l >> RNA_counts
#echo DONE $(date)

echo STARTED $(date)
touch DNA_counts
echo "DNA_1" >> DNA_counts
find ./ -iname \*DNA_1\* -print0 | xargs -0 -I file egrep '^@.*' file \
   | wc -l >> DNA_counts
echo "Counting reads in DNA_2"
echo "DNA_2" >> DNA_counts
find ./ -iname \*DNA_2\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> DNA_counts
echo "Counting reads in DNA_3"
echo "DNA_3" >> DNA_counts
find ./ -iname \*DNA_3\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> DNA_counts
echo "Counting reads in DNA_4"
echo "DNA_4" >> DNA_counts
find ./ -iname \*DNA_4\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> DNA_counts
echo DONE $(date)
