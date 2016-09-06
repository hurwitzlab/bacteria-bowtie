#!/usr/bin/env bash

if [[ -e ./mini-config-mouse.sh ]]; then
    source ./mini-config-mouse.sh
else
    echo "no config-mouse"
    exit 1
fi

echo "Making id_to_gene.tab"
sleep 1

perl -nle '($id)=/gene_id "([^"]*)/; ($gene)=/gene_name "([^"]*)/; print "$id\t$gene"' \
    $RWORK_DIR/genes.gtf > $RWORK_DIR/id_to_gene.tab
    
echo "Leaving directory"
