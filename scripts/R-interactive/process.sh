#!/usr/bin/env bash

if [[ -e ./mini-config.sh ]]; then
    source ./mini-config.sh
else
    echo "no config"
    exit 1
fi

echo "Munging gff file"
sleep 5

perl -nle '($id)=/ID=([^;]*)/; ($prod)=/product=([^;]*)/; print "$id\t$prod"' \
    $RWORK_DIR/all_fixed.gff > $RWORK_DIR/id_to_product.tab

echo "Leaving directory"
