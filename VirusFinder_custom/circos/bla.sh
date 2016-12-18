awk '{ gsub ("chr","hs", $1); gsub ("\"","", $18); gsub (";","", $18); print $1,$4,$5,$18 }' /media/StorageOne/HTS/PublicData/hg18_Human_genes.gtf > ./data/text.genes.hg18.txt
sed -i '/^#/d' ./data/text.genes.hg18.txt
