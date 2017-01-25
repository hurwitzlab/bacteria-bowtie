#!/usr/bin/env bash

if [[ -e ./mini-config-combined.sh ]]; then
    source ./mini-config-combined.sh
else
    echo "no config"
    exit 1
fi

echo "Making id_to_product.tab [description]"
sleep 1

sed -ibak s/Name=/product=/g $RWORK_DIR/combined.gff

perl -nle '($id)=/ID=([^;]*)/; ($prod)=/product=([^;]*)/; print "$id\t$prod"' \
    $RWORK_DIR/combined.gff > $RWORK_DIR/id_to_product.tab

echo "Making id_to_gene.tab"
sleep 1

perl -nle '($id)=/ID=([^;]*)/; ($gene)=/gene=([^;]*)/; print "$id\t$gene"' \
    $RWORK_DIR/combined.gff > $RWORK_DIR/id_to_gene.tab

echo "Leaving directory"
