#! /bin/bash

## update 30.11.2023
snpFile=$1
key=$2
interactions=$3

seed=1 # for sneep run
rounds=0

## get the enhancer from dennis and remove those interactions not linked to a gene
#python parseAbcInteractions.py ${interactions} ABC_interaction_updated_19122023.txt ## is done in workflow_3.0.sh file 

## seperated them based on their celltype
python seperateInteractions.py ABC_interaction_updated_22052024.txt  ABC_interactions_seperated_updated_22052024/

for filename in  ABC_interactions_seperated_updated_22052024/*; do

	#echo ${filename}
	IFS='_' read -r -a array <<< "$filename" # split string into array, sep = _
	helper=${array[6]} 
	size=${#helper}-4
	celltype=${helper:0:${size}}
	echo "$celltype"

	echo ${seed}

	time src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_${celltype}_${key}/ -n 20  -p 0.5 -c 0.001 -b frequence.txt   -f ${filename} -r ${filename}  -g geneId_geneName.txt   -j 0 -l ${seed} -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt ${snpFile} hg38.fa  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt

	rm -r sneep_cadlincGWAS_${celltype}_${key}/PFMs/
	seed=$((seed+1))
	echo ${seed}

done
