#!/usr/bin/env bash

if [[ -e ./mini-config-unknown.sh ]]; then
    source ./mini-config-unknown.sh
else
    echo "no config"
    exit 1
fi

echo "Making id_to_product.tab [description]"
sleep 1

sed -ibak s/Name=/product=/g $RWORK_DIR/unknown.gff

perl -nle '($id)=/ID=([^;]*)/; ($prod)=/product=([^;]*)/; print "$id\t$prod"' \
    $RWORK_DIR/unknown.gff > $RWORK_DIR/id_to_product.tab

echo "Leaving directory"
