#!/bin/sh

infasta=$1
outfasta=$2

cat $infasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ print $1"@len="length($2)"@path=[0:0-000@000:000-000@0000:000-000]",$2}' > $outfasta
sed -i 's/ /\n/g' $outfasta
sed -i 's/@/ /g' $outfasta

