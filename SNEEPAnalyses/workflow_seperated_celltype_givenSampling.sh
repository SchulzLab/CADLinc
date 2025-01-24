#!/usr/bin/env bash
snpFile=$1
key=$2
background=$3

seed=1 # for sneep run
rounds=500

output="sneep_result_celltypeSpecific_withBackground/countRegREMsPerCelltype/"

mkdir ${output}

for filename in  ABC_interactions_seperated_updated_19122023/*; do

	IFS='_' read -r -a array <<< "$filename" # split string into array, sep = _
	helper=${array[6]} 
	size=${#helper}-4
	celltype=${helper:0:${size}}
	echo "$celltype"

	time differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_${celltype}_${key}/ -n 20  -p 0.5 -c 0.001 -b frequence.txt -f ${filename} -r ${filename}  -g geneId_geneName.txt -i ${background}/randomSNPs_${celltype}/randomSNPs_sneep/  -j ${rounds} -l ${seed} -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt ${snpFile} hg38.fa allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt

done
