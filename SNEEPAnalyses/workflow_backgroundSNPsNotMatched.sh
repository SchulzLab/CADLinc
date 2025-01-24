#! /bin/bash

## update 30.11.2023
snpFile=$1
key=$2
interactions=$3

seed=1 # for sneep run
rounds=1

for filename in  ABC_interactions_seperated_updated_19122023/*; do

	#echo ${filename}
	IFS='_' read -r -a array <<< "$filename" # split string into array, sep = _
	helper=${array[6]} 
	size=${#helper}-4
	celltype=${helper:0:${size}}
	echo "$celltype"

	echo ${seed}

	echo "	time differentialBindingAffinity_multipleSNPs -o randomSNPs_cadlincGWAS_2.52x0.00001_22_02/randomSNPs_${celltype}/sneep_cadlincGWAS_${celltype}_${key}/ -n 20  -p 0.5 -c 0.001 -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt  -f ${filename} -r /projects/cadlinc/work/SNEEP_analysis/${filename}  -g geneId_geneName.txt   -j ${rounds} -l ${seed} -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt ${snpFile}/randomSNPs_${celltype}/notMatchedSNPs_${celltype}.txt hg38.fa"

	seed=$((seed+1))
	echo ${seed}

done
