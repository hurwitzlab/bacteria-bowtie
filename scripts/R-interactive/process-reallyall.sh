#!/usr/bin/env bash

if [[ -e ./mini-config-reallyall.sh ]]; then
    source ./mini-config-reallyall.sh
else
    echo "no config"
    exit 1
fi

echo "Making id_to_product.tab [description]"
sleep 1

perl -nle '($id)=/ID=([^;]*)/; ($prod)=/product=([^;]*)/; print "$id\t$prod"' \
    $RWORK_DIR/all-refseq-CDS.gff > $RWORK_DIR/id_to_product.tab

echo "Making id_to_gene.tab"
sleep 1

perl -nle '($id)=/ID=([^;]*)/; ($gene)=/gene=([^;]*)/; print "$id\t$gene"' \
    $RWORK_DIR/all-refseq-CDS.gff > $RWORK_DIR/id_to_gene.tab

echo "Removing extra headers"
cd $RWORK_DIR
head -1 total.fpkm_tracking > header
egrep -v '^tracking_id.*' total.fpkm_tracking > temp01
cat header temp01 > total.fpkm_tracking
rm header
rm temp01

echo "Leaving directory"
