export IFS=$'\n'

for SP in $(cat ./LPS_search_list); do
    if [[ $SP =~ "^#" ]]; then
        continue
    else
        echo Doing "$SP"
        egrep "$SP" all-refseq-CDS.gff >> LPS-refseq-CDS.gff
    fi
done

export IFS=$' \t\n'
