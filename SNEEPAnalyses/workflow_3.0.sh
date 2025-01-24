#! /bin/bash

## Step1: extract significant SNPs (p-value <= 0.000001) from FINAL_1MH_CAD_GWAS_summary_stats.tsv   
## important to use the summary statistic with the already swapped alleles (more details are given in workflow_allSNPsSummary.sh line 5)
#awk '{if($9 <= 2.52*0.00001) print  $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t-\t" $9}' GWAS_summary_NOT_Filtered_betaSwap_hg19.txt  > GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_09_02.txt

## add header #chr    start   end     allele1 allele2 rsID    pvalueSummaryStatistic(column:Pvalue) to GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_09_02.txt
## get proxy SNPs from 100 genome prject
python getProxySNPs_2.0.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_09_02.txt  0.8 snipa_data/ GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_16_02_lead_proxySNPs.txt

## check direction of the proxy SNPs
python swapProxySNPs.py FINAL_1MH_CAD_GWAS_summary_stats.tsv  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_16_02_lead_proxySNPs.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_lead_proxySNPs_swapped.txt 

##changed so that also the leadUnknown snps have an id XX123
python getUniqEntries_hg19.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_lead_proxySNPs_swapped.txt GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_lead_proxySNPs_swapped_uniq.txt

## with adapted lifting
python liftPositions_tohg38_2.0.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_lead_proxySNPs_swapped_uniq.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt

## get all LD SNPs
## sort file (starting from 4 letter in first column (ignore chr)
sort -k1.4n GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt  > GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sorted.txt

python getProxySNPs_proxy.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sorted.txt 0.8 snipa_data/ GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_allLDSNVs.txt

cut -f1,2,3,4,5,6,7  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt >  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt


## intersection only the lead SNPs with ABC interactions 
grep "lead"  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt > GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_ONLY_lead_swapped_uniq.txt

## intersection all GWAS SNP with ABC interactions
bedtools intersect -a GWAS_summary_NOT_Filtered_betaSwap_hg38.txt -b ABC_interaction_updated_19122023_merged.txt -u > intersection_GWAS_summary_NOT_filtered_ABC_interactions.bed

## rerun with update enhancer-gene links
python parseAbcInteractions.py  enhancer_interactions_avghg38_unique/ ABC_interaction_updated_22052024.txt
## newest version with sca;e file as required argument
time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_22_05/ -n 20  -p 0.5 -c 0.001 -b frequence.txt   -r ABC_interaction_updated_22052024.txt      -g geneId_geneName.txt  -j 0 -l 23 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt

## with uniq peak files -> keep in mind to use the seperatedPFMs_CADversion.py file instead of the newest one, we used the sneep evrsion where the PWMs were computed with the varying epsilon
bash workflow_seperated_celltype.sh GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt  2.52x0.00001_noSampling_27_05_2024 ../enhancer_interactions_avghg38_unique/

## call snipsnap to get the new background snps
## get snps in snpsnap format
python getSNISNAP_input.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_16_02_lead_proxySNPs_swapped_uniq.txt notMappedSNPs_2.52x0.00001_22_02.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_snpsnap.txt

## call snpsnap 
bash run_snpsnap.sh  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_22_02_snpsnap.txt  GWAS_cadlinc_2.52x0.00001_hg19_22_02

## get the matchedSNPs per celltype and those which are not matched and for which we need to run sneep
python parseSNPSNAP_output.py GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt ../enhancer_interactions_avghg38/celltypes.txt  snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/matched_snps.tsv  snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/match_quality_per_snp.tsv randomSNPs_cadlincGWAS_2.52x0.00001_22_02/

##run sneep for the not matched snps with background sampling
bash workflow_backgroundSNPsNotMatched.sh randomSNPs_cadlincGWAS_2.52x0.00001_22_02/  2.52x0.00001_snpSampling_22_02 ../enhancer_interactions_avghg38/ >output_22_02_2.52_randomSNPs_sneep.txt 

#get allele info of the matched snps and lift them to hg38
## parse snipa file to only keep the information per snp once 
time python parseSNIPA_SNPS.py snipa_data/ snipa_data/snipa_snps_with_allele_info.txt

## get alleles for the snps
python getRandomSNPsSNEEP.py snipa_data/snipa_snps_with_allele_info_sorted.txt randomSNPs_cadlincGWAS_2.52x0.00001_22_02/

## lift matched snpsnap snps to hg38, split them into the backgroudn rounds and add the sampled sneep snps 
## redo for pDC and mDC all other celltypes are the same as before -> data for pDC and mDC is copied to the same folder than before
python splitRandomSNPs.py randomSNPs_cadlincGWAS_2.52x0.00001_22_02/ 

## run sneep with background snpsnap snps
## rerun because of the updated sneep version 
time bash workflow_seperated_celltype_givenSampling.sh GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt 2.52x0.00001_backgroundSampling_22_02 randomSNPs_cadlincGWAS_2.52x0.00001_22_02/ > output_22_02_seperatedCellTypes_withBackground.txt

## plot the data 

##derive odds-ratio
##adapted pseudocount
python3 computeOddsRatio.py  sneep_result_celltypeSpecific_withBackground/ TF_names.txt  numberSNPsPerCelltype.txt heatmapInput_oddsRatio_cutoffs_pseudocount.txt heatmapInput_lossGain_oddsRatio_cuto_pseudcount.txt

## get the background snps for running sneep with mergen enhnacer-gene interactions over all celltypes
# need to modify parseSNPSNAP_output.py -> parseSNPSNAP_output_mergedEnhancers.py
mkdir randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/

bedtools merge -i <(sort -k1,1 -k2,2n ABC_interaction_updated_19122023.txt ) > ABC_interaction_updated_19122023_merged.txt  

time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_27_02/ -n 20  -p 0.5 -c 0.001 -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt -f ABC_interaction_updated_19122023_merged.txt    -r ABC_interaction_updated_19122023.txt      -g geneId_geneName.txt  -j 0 -l 23 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa

python parseSNPSNAP_output_mergedEnhancers.py  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt  sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_27_02/ snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/matched_snps.tsv  snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/match_quality_per_snp.tsv randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/

##rerun sneep to get the random snps for notMatchedSNP of snpsnap
time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/randomSNPs/sneep_cadlincGWAS_2.52x0.00001_snpSampling_unifiedCelltypes_27_02/ -n 20  -p 0.5 -c 0.001  -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt   -g geneId_geneName.txt   -j 500 -l 123 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/randomSNPs/notMatchedSNPs.txt hg38.fa

python getRandomSNPsSNEEP.py snipa_data/snipa_snps_with_allele_info_sorted.txt randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/

# split in 500 background files
python splitRandomSNPs_unifiedCelltypes.py randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/ sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_27_02/

## run sneep with background analaysis
time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_backgroundSampling_27_02/ -n 20  -p 0.5 -c 0.001 -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt  -f ABC_interaction_updated_19122023_merged.txt -r ABC_interaction_updated_19122023.txt -g geneId_geneName.txt -i randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02/randomSNPs/randomSNPs_sneep/  -j 500 -l 12 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa

## run sneep with background analysis and merged REMs 
## get uniq rem-gene interactions per gene over all celltypes
bash workflow_uniqRemGeneInteractions.sh ABC_interaction_updated_19122023.txt ABC_interactionsMergedPerGene/ > mergeREMs.log

time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_backgroundSampling_10_05/  -n 20  -p 0.5 -c 0.001 -b frequence.txt   -f ABC_interaction_updated_19122023_merged.txt -r ABC_interactionsMergedPerGene/ABC_interactionsMergedGenes.txt -g geneId_geneName.txt -i randomSNPs_cadlincGWAS_2.52x0.00001_unifiedCelltypes_27_02_mergedEnhancer/randomSNPs/randomSNPs_sneep/  -j 500 -l 12 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt &>sneep.log

## run sneep with background sampling and epiregio rems

## add pseudocount from 1
python computeOddsRatio_mergedEnhancer.py sneep_cadlincGWAS_2.52x0.00001_mergedEnhancers_backgroundSampling_27_02/ TF_names.txt  10880  oddsRatio_mergedEnhancer.txt lossGain_mergedEnhancer.txt

## rerun for the updated pDC and mDC peak file
##instead of oddsRatio_mergedEnhancer.txt use heatmapInput_oddsRatio_cutoffs_pseudocount_labeled.txt and sneep_cadlincGWAS_2.52x0.00001_22_05/result for target genes to get the updated once for pDC and mDC
python getTF_evidence.py TF_table_2905.txt  TF_evidence/  JASPAR2022/TF_ensemblID_name_human_JASPAR2022_GRCh38p13.txt  FisherTestCoExpressedGenes_sneep_geneBased_100_greater_JointKnownGenes_2509.txt 2509_JointKnownGenes_GeneTable.txt  TF_evidence/overlapCAD_genes_0810.txt TF_evidence/overlapingleCellData_0810.txt TF_table_0810_extended.txt

## with fixed pDC and mDC peak files and newest gene list from zhifen
python getTF_evidence.py TF_table_2905.txt  TF_evidence/  JASPAR2022/TF_ensemblID_name_human_JASPAR2022_GRCh38p13.txt  FisherTestCoExpressedGenes_sneep_geneBased_100_greater_fixedDC_0409.txt ../230524_FixDC_GeneTable.txt  TF_evidence/overlapCAD_genes_0409.txt TF_evidence/overlapingleCellData_0409.txt TF_table_0409_extended.txt

## check which eQTLs are linked to the same genes than our prediction
python check_eQTL_data.py  sneep_result_celltypeSpecific/ eQTL_TF_analysis/All_eQTL_GTEx_SNEEP_2.52e05_TF_SNP_Celltype.csv  eQTL_TF_analysis/eQTL_STARNET_SNEEP_2.52e05_TF_SNP_Celltype.txt OpenGenesCell.pkl SNPs_Enhancer_links_eQTLs.txt

##
## repeat  for all 1%FDR cutoff snps without enhacer gene links
time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_noABC_Interactions_27_02/ -n 20  -p 0.5 -c 0.001 -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa

mkdir randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/

python parseSNPSNAP_output_mergedEnhancers.py  GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq.txt sneep_cadlincGWAS_2.52x0.00001_noABC_Interactions_27_02/ snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/matched_snps.tsv  snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_22_02_matched_snps_for_input_snps/match_quality_per_snp.tsv randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/

time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/randomSNPs/sneep_cadlincGWAS_2.52x0.00001_snpSampling_noABC_interactions_27_02/ -n 20  -p 0.5 -c 0.001  -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt   -g geneId_geneName.txt   -j 500 -l 123 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/randomSNPs/notMatchedSNPs.txt hg38.fa

python getRandomSNPsSNEEP.py snipa_data/snipa_snps_with_allele_info_sorted.txt randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/

python splitRandomSNPs_unifiedCelltypes.py randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/ sneep_cadlincGWAS_2.52x0.00001_noABC_Interactions_27_02/

time pipelineSNEEP/src/differentialBindingAffinity_multipleSNPs -o sneep_cadlincGWAS_2.52x0.00001_noABC_Interactions_backgroundSampling_27_02/  -n 20  -p 0.5 -c 0.001 -b frequence.txt  -s  allele_specific_SNPs/scalesPerMotif/estimatedScalesPerMotif_1.9.txt  -i randomSNPs_cadlincGWAS_2.52x0.00001_noABC_interactions_27_02/randomSNPs/randomSNPs_sneep/  -j 500 -l 12 -k dbSNP/dbSNP_2022_11_16/dbSNPs_sorted.txt  combined_Jaspar_Hocomoco_Kellis_human_transfac_jaspar2022_withoutCTCF_MA12929.txt GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_sneep.txt hg38.fa


python computeOddsRatio_mergedEnhancer.py sneep_cadlincGWAS_2.52x0.00001_noABC_Interactions_backgroundSampling_27_02/ TF_names.txt  72433  oddsRatio_noABC_interactions.txt lossGain_noABC_interactions.txt

bash workflow_seperated_celltype.sh GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_16_02_lead_proxySNPs_swapped_uniq_sneep.txt  2.52x0.00001_noSampling_16_02 ../enhancer_interactions_avghg38/


