#!/bin/bash

SNPsPerGene=$1
outputDir=$2

## Calculate pair-wise ld correlation for 1000 genome per cromosome
## Results are save in folder LD


for i in $(ls alkesgroup.broadinstitute.org/LDSCORE/GRCh38/plink_files/*.fam);
do
	(filename=$(basename $i)
	file="${filename%%.fam}"
	output_file=${file#"1000G.EUR.hg38."}
	echo $filename
	echo $file
	echo $output_file
	echo "--------------"

	## for gene in geneFile
	for gene in $(ls ${SNPsPerGene}/${output_file}/*.txt);

	do

		gene_name=$(basename $gene)
		## read all SNPs for this Gene
		snp_string=$(paste -sd, $gene)
		echo $snp_string

		/projects/cadlinc/work/gene_based_test/plink-1.07-x86_64/plink --noweb --bfile /projects/cadlinc/work/gene_based_test/alkesgroup.broadinstitute.org/LDSCORE/GRCh38/plink_files/${file} --snps $snp_string --r2 --matrix --ld-snp-list ${gene} --out ${outputDir}/${gene_name}

	done) & ## runs in parallele mode defined by ()& and wait

done

wait

