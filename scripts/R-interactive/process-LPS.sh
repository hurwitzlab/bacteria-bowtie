#!/usr/bin/env bash

if [[ -e ./mini-config-LPS.sh ]]; then
    source ./mini-config-LPS.sh
else
    echo "no config"
    exit 1
fi

echo "Making id_to_product.tab [description]"
sleep 1

perl -nle '($id)=/ID=([^;]*)/; ($prod)=/product=([^;]*)/; print "$id\t$prod"' \
    $RWORK_DIR/LPS-refseq-CDS.gff > $RWORK_DIR/id_to_product.tab

echo "Making id_to_gene.tab"
sleep 1

perl -nle '($id)=/ID=([^;]*)/; ($gene)=/gene=([^;]*)/; print "$id\t$gene"' \
    $RWORK_DIR/LPS-refseq-CDS.gff > $RWORK_DIR/id_to_gene.tab
    
echo "Leaving directory"
