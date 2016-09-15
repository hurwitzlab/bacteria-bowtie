echo STARTED $(date)
touch RNA_counts
echo "Counting reads in RNA_1"
echo "RNA_1" >> RNA_counts
find ./ -iname \*RNA_1\* -print0 | xargs -0 -I file egrep '^@.*' file \
   | wc -l >> RNA_counts
echo "Counting reads in RNA_2"
echo "RNA_2" >> RNA_counts
find ./ -iname \*RNA_2\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> RNA_counts
echo "Counting reads in RNA_3"
echo "RNA_3" >> RNA_counts
find ./ -iname \*RNA_3\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> RNA_counts
echo "Counting reads in RNA_4"
echo "RNA_4" >> RNA_counts
find ./ -iname \*RNA_4\* -print0 | xargs -0 -I file egrep '^@.*' file \
    | wc -l >> RNA_counts
echo DONE $(date)
