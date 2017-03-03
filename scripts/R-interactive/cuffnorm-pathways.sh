#!/usr/bin/env bash

if [[ -e ./mini-config-cuffnorm.sh ]]; then
    source ./mini-config-cuffnorm.sh
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

#We need the ec_number
#So we can fetch the pathways using keggFind
#e.g. keggFind("ko","[EC:1.2.1.3]")->thing
#and keggLink("pathway",names(thing))->thing2

echo "Making id_to_ecnumber.tab"
sleep 1
perl -nle '($id)=/ID=([^;]*)/; ($ec)=/ec_number=([^;]*)/; print "$id\t$ec"' \
    $RWORK_DIR/all-refseq-CDS.gff > $RWORK_DIR/id_to_ecnumber.tab

echo "Leaving directory"
