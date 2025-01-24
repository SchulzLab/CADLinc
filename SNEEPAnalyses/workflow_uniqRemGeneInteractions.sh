#!/usr/bin/env bash

remGeneInteractions=$1 #ABC_interaction_updated_19122023.txt  over akk celltypes not merged
outputDir=$2

## get genes of interest

while read g;
do
	echo "${g}"
	echo "grep -P "${g}\t" ${remGeneInteractions} > ${outputDir}/rems_${g}.txt"
	grep -P "${g}\t" ${remGeneInteractions} > ${outputDir}/rems_${g}.txt

	echo "sort -k1,1 -k2,2n ${outputDir}/rems_${g}.txt > ${outputDir}/rems_sorted_${g}.txt"
	sort -k1,1 -k2,2n ${outputDir}/rems_${g}.txt > ${outputDir}/rems_sorted_${g}.txt

	echo "bedtools merge -i ${outputDir}/rems_sorted_${g}.txt > ${outputDir}/rems_merged_${g}.txt"
	bedtools merge -i ${outputDir}/rems_sorted_${g}.txt > ${outputDir}/rems_merged_${g}.txt
	echo "awk -v var="$g" '{print $1 "\t" $2 "\t" $3 "\t" var "\tallCelltypes\t-\t-\t-\t-\t-\t-\t-"}' ${outputDir}/rems_merged_${g}.txt > ${outputDir}/rems_merged_parsed_${g}.txt"
	awk -v var="$g" '{print $1 "\t" $2 "\t" $3 "\t" var "\tallCelltypes\t-\t-\t-\t-\t-\t-\t-"}' ${outputDir}/rems_merged_${g}.txt > ${outputDir}/rems_merged_parsed_${g}.txt

	rm ${outputDir}/rems_merged_${g}.txt
	rm ${outputDir}/rems_sorted_${g}.txt
	rm ${outputDir}/rems_${g}.txt

done < ${outputDir}/genes.txt

## cat over all 60.000 resulting files
## run in output file
time find . -type f -name "rems_merged_parsed_ENSG00000*" -print0 | xargs -0 cat > ABC_interactionsMergedGenes.txt



