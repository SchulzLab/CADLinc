#!/bin/bash

input_snps_file=$1 # "GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_16_02_snpsnap.txt"
analysis_name=$2 #"GWAS_cadlinc_2.52x0.00001_hg19_16_02"

## parameters that are for all runs the same
n_matches=800
ldbud_r2="friends_ld05"
output_root="snpsnap_results/"
db_file="snpsnap/v2_2020_04_27/ld0.8_collection.tab.gz" ## the collection is from qianqian and she downloaded it from the snpsnap webpage when it was still working (currently it is not working)

## actiavte conda enviroment
#conda activate SNIPSNAP  
## create the output dir
mkdir ${output_root}

python match_snps.py $input_snps_file $n_matches $ldbud_r2 $db_file $analysis_name $output_root

