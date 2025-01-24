#!/bin/bash

### download of the PLINK LD data orginal from here https://alkesgroup.broadinstitute.org/LDSCORE/GRCh38/ (identified based on nikolettas file structure)
### now moved to https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays/LDSCORE/GRCh38?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false

## : cat ../../alkesgroup.broadinstitute.org/LDSCORE/GRCh38/plink_files/1000G.EUR.hg38.*.bim > PLINK_SNPs.txt

############
### STEP 1: lift SNV from the GWAS summary statistic, which is hg19, to hg38
### takes as input summary statistics and return lifted position to hg38, header and additional info is the same as before
############

python liftHG19_toHG38.py FINAL_1MH_CAD_GWAS_summary_stats.tsv GWAS_SUMMARY_hg38.txt
## lost snps when lifting to hg38: 7060

############
### STEP 2: get ABC interactions 
############

### peak location will also be removed from the script
### only keep chr start end gene 0-based 

cut -f1,2,3,4 ABC_interaction_updated_19122023.txt >  ABC_interactions_updated_19122023_genes.txt

############
### STEP 3: intersect GWAS summary SNVs with gene body from annotation file 
### input file for 
############

### write 0-based SNV coordinates as columns (to get bed format) followed by info from the GWAS file without header
awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print "chr"$2,$3-1,$3,$0 }}' GWAS_SUMMARY_hg38.txt  > GWAS_SUMMARY_hg38_positions.txt

### use bedtools intersect to get the overlap to the annotation file, output files holds info annotation file followed by 0-based coordinates, GWAS info and overlap (1)
bedtools intersect -wo -a gencode.v38.annotation.gtf -b GWAS_SUMMARY_hg38_positions.txt > SNPs_overlapping_annotation.txt

### filter the overlap for gene and extract relevant information 
### pos 9 contains the information gene_id "ENSG00000240361.2"; transcript_id "ENST00000642116.1"; gene_type ... extract only gene_id 
### keep all information from GWAS file and add gene id as last column
awk '$3=="gene"' SNPs_overlapping_annotation.txt | awk 'BEGIN{FS="\t";OFS="\t"}{split($9,a,";")};{split(a[1],b,"\"")};{split(b[2],c,".")};{ print $13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,c[1]}' > SNPstoGenes.txt

############
### STEP 4: intersect GWAS summary SNVs with REMs
############

bedtools intersect -a GWAS_SUMMARY_hg38_positions.txt -b ABC_interactions_updated_19122023_genes.txt  -wo > SNPs_overlapping_REMs.txt
### results in 0-based position, GWAS info and	chr1	804896	805046	ENSG0000023009 1
### remove 0-based position, keep GWAS info and add at last column the gene id
awk 'BEGIN{FS="\t";OFS="\t"}{ print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$28}' SNPs_overlapping_REMs.txt > SNPstoREMs.txt

############
#### run gene-based test only with rSNVs (from all GWAS SNVs) in REMs
############

## with sneep I already computed rSNVs without a filter chunkedGWAS_summary_NOT_Filtered_betaSwap_hg38_sneep/positionAllRegSNPs_sorted.txt -> not overlapped with REMs
## bring in bed-like format and overlap with REMs (using vim replace : by tab und - by tab) -> positionAllRegSNPs_sorted.bed

## get info from GWAS_SUMMARY for each SNP
bedtools intersect -wo -a GWAS_SUMMARY_hg38_positions.txt -b chunkedGWAS_summary_NOT_Filtered_betaSwap_hg38_sneep/positionAllRegSNPs_sorted.bed> regSNVs_overlapping_REMs.txt  ## 69.52.103 number TF_SNVs

## get overlap with REMs
bedtools intersect -wo -a regSNVs_overlapping_REMs.txt -b ABC_interactions_updated_19122023_genes.txt  > regSNVs_overlappingREMs_.txt ## 22.143.936 number TF-SNVs overlapping a REM

awk 'BEGIN{FS="\t";OFS="\t"}{ print $4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$32}' regSNVs_overlappingREMs_.txt > regSNVstoREMs.txt

## keep all genebody SNVs and only add rSNVs  -> keep_regSNVs_REMs.py
############
### STEP 5: prepare the SNV sets (gene-body, enhancers, gene-body with enhancer) for gene-based test 
############

##header is from GWAS file plus last column GENE
# add header and uniq
# only SNVs in genebody: SNPstoGenes_uniq.txt
(echo -e "MarkerName\tCHR\tBP\tAllele1\tAllele2\tFreq1\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tCases\tEffective_Cases\tN\tMeta_analysis\tGENE"; sort SNPstoGenes.txt | uniq )  > SNPstoGenes_uniq.txt

############
### STEP 6: add rsID from PLINK data 
############

## PLINK_SNPs are hg38 (checked via dbSNP database for a few entries)  PLINK_SNPs.txt

## for SNVs in gene body
python assign_names_SNPs.py PLINK_SNPs.txt SNPstoGenes_uniq.txt SNPstoGenes_rsID.txt # #SNVs 6.655.529

############
### STEP 7: seperate rsIDs per gene
############

## sort based on the gene is necessary for seperating the rsID per genes 
sort -k28 SNPstoGenes_rsID.txt  > SNPstoGenes_rsID_sorted.txt 

mkdir SNPsperGene_genebody/
$END=22
for i in $(seq 1 $END); do
	echo $i;
	mkdir SNPsperGene_genebody/${i}/
done

## seperate rsIDs per gene
python group_SNP_genes.py  SNPstoGenes_rsID_sorted.txt  SNPsperGene_genebody/ genes_genebody.txt  ## output: number genes: 48187 and length list of genes: 48187, counter SNPs considered: 6649019, counter SNPs duplicated rsID but differ in alleles: 6510

## check with find . -type f | wc -l if the number of gene files is identical to the number of genes (checked and is fine)

############
### STEP 8: get LD structure
############

mkdir LD_genebody
time bash LD_calculation.sh SNPsperGene_genebody/  LD_genebody/ &> log_LD_genebody.log

############
### STEP 9: run gene-based test
############

# for all genes for the genebody
time Rscript GATES.R  SNPstoGenes_rsID.txt genes_genebody.txt LD_genebody/ result_gene_based_test_genebody.txt &> gene_based_test_genebody.log ## computes gene based test for 43.137 genes (from 48.187 genes)
python parseOutputGeneBasedTest.py  result_gene_based_test_genebody.txt gencode.v38.annotation.gtf  result_gene_based_test_genebody_sig.txt


