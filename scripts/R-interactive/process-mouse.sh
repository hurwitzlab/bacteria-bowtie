#!/usr/bin/env bash

if [[ -e ./mini-config-mouse.sh ]]; then
    source ./mini-config-mouse.sh
else
    echo "no config-mouse"
    exit 1
fi
#Y          RefSeq  exon    2397909 2397997 .       +       .       gene_id "Gm3376"; gene_name "Gm3376"; p_id "P8059"; transcript_id "XM_001475806.2"; tss_id "TSS18259";

echo "Making id_to_gene.tab"
sleep 1

perl -nle '($id)=/transcript_id "([^"]*)/; ($gene)=/gene_name "([^"]*)/; print "$id\t$gene"' \
    $RWORK_DIR/mouse-exons.gtf > $RWORK_DIR/id_to_gene.tab
    
echo "Leaving directory"
