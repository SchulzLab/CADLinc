library(dplyr)
library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(extrafont)
library(scales)

# Process coloc output files
# Data preparation
GTEx_coloc_significant_SNPs_08 <- read.csv("results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_significant_results_PPH4_08_Pval_5e_8.csv")
GTEx_coloc_significant_SNPs_08$Study <- 'GTEx'

STARNET_coloc_significant_SNPs_08 <- read.csv("results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_STARNET_significant_results_PPH4_08_Pval_5e_8.csv")
STARNET_coloc_significant_SNPs_08$Study <- 'STARNET'

# Merge 2 coloc outputs and select only similar columns
# Make a subset
coloc_output_GTEx_08 <- GTEx_coloc_significant_SNPs_08[, c("SNP", "position", "SNP.PP.H4", "gene_ensembl", "chr", "position_hg19", "chr_bp_position_hg19", "slope", "p_value_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", "SE", "sample_size", "tissue", "Trait", "test", "Gene.name",  "Biotype", "Sample_size_eQTL", "Study")]
coloc_output_STARNET_08 <- STARNET_coloc_significant_SNPs_08[, c("SNP", "position", "SNP.PP.H4", "gene_ensembl", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "p_value_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", "SE", "sample_size", "tissue", "Trait", "test", "Gene.name",  "Biotype", "Sample_size_eQTL", "Study")]

# Rename some columns
names(coloc_output_GTEx_08)[names(coloc_output_GTEx_08) == "slope"] <- "beta_eQTL"
coloc_output_08 <- rbind(coloc_output_GTEx_08, coloc_output_STARNET_08)

# Add new column with 2 types in biotype
coloc_output_08 <- coloc_output_08 %>%
  mutate(Biotype_binary = ifelse(Biotype == "protein_coding", "protein-coding", "non-coding"))
table(coloc_output_08$Biotype_binary)

# Rename traits for plotting
coloc_output_08$Trait_labels <- coloc_output_08$Trait
coloc_output_08$Trait_labels <- revalue(coloc_output_08$Trait_labels,
                                   c("Atrial_fibrillation" = "Atrial fibrillation",
                                     "Heart_Failure" = "Heart failure",
                                     'Resting_heart_rate' = 'Resting heart rate',
                                     "Type_2_diabetes" = "Type 2 diabetes",
                                     "NAFLD" = "Nonalcoholic fatty liver disease",
                                     "HDL_cholesterol" = "HDL cholesterol",
                                     "LDL_cholesterol" = "LDL cholesterol",
                                     "Triglycerides" = "Triglycerides",
                                     "Total_cholesterol" = "Total cholesterol",
                                     "Non-HDL_cholesterol" = "Non-HDL cholesterol",
                                     "Fasting_glucose" = "Fasting glucose",
                                     "Fasting_insulin" = "Fasting insulin",
                                     "HbA1c" = "HbA1c",
                                     "Diastolic_blood_pressure" = "Diastolic blood pressure",
                                     "Systolic_blood_pressure" = "Systolic blood pressure",
                                     "Pulse_pressure" = "Pulse pressure",
                                     "Two_hour_glucose" = "Two hour glucose test",
                                     "Basophil_count" = "Basophil count",
                                     "Eosinophil_count" = "Eosinophil count",
                                     "Lymphocyte_count" = "Lymphocyte count",
                                     "Monocyte_count" = "Monocyte count",
                                     "Neutrophil_count" = "Neutrophil count",
                                     "Red_blood_cell_count" = "Red blood cell count",
                                     "Platelet_count" = "Platelet count",
                                     "White_blood_cell_count" = "White blood cell count",
                                     "Body_mass_index" = "Body mass index",
                                     "Waist-hip_ratio" = "Waist-hip ratio",
                                     "pad_primary" = "Peripheral arterial disease",
                                     "pad_diabetes" = "Peripheral arterial disease",
                                     "pad_nodiabetes" = "Peripheral arterial disease",
                                     "pad_eversmoker" = "Peripheral arterial disease",
                                     "pad_neversmoker" = "Peripheral arterial disease"))

head(coloc_output_08)

# Rename Tissue for plotting, merge them in groups
coloc_output_08$tissue_labels <- coloc_output_08$tissue
coloc_output_08$tissue_labels <- revalue(coloc_output_08$tissue_labels,
                                   c('Adipose_Subcutaneous' = 'Adipose tissue',
                                     "Adipose_Visceral_Omentum" = "Adipose tissue",
                                     "SF" = "Adipose tissue",
                                     "VAF" = "Adipose tissue",
                                     'AOR' = 'Artery',
                                     "Artery_Aorta" = "Artery",
                                     "Artery_Coronary" = "Artery",
                                     "Artery_Tibial" = "Artery",
                                     "MAM" = "Artery",
                                     "BLD" = "Blood",
                                     "Whole_Blood" = "Blood",
                                     "Heart_Atrial_Appendage" = "Heart",
                                     "Heart_Left_Ventricle" = "Heart",
                                     "Kidney_Cortex" = "Kidney cortex",
                                     "LIV" = "Liver",
                                     "Liver" = "Liver",
                                     "Muscle_Skeletal" = "Skeletal muscle",
                                     "SKLM" = "Skeletal muscle",
                                     "Nerve_Tibial" = "Tibial nerve",
                                     "Spleen" = "Spleen"
                                     ))
head(coloc_output_08)

# Add the information regarding CAD known/novel candidate genes
gene_table <- read.csv('2509_JointKnownGenes_GeneTable.txt', sep = '\t', header = TRUE)
gene_table$known_CAD_loci_gene_new_list <- ifelse(
  gene_table$known.CAD.loci.gene == 'True' | 
  gene_table$Coloc.GTEx == 'True' | 
  gene_table$Coloc.STARNET == 'True', 
  "True", 
  "False"
)

# Add the columns to the coloc gene table, perform a left join
# 'Conserved.in.mouse' column
coloc_output_08 <- coloc_output_08 %>%
  left_join(gene_table %>% dplyr::select(Gene.name, Conserved.in.mouse), by = "Gene.name")

# 'known.CAD.loci.gene' column
coloc_output_08 <- coloc_output_08 %>%
  left_join(gene_table %>% dplyr::select(Gene.name, known_CAD_loci_gene_new_list), by = "Gene.name")
head(coloc_output_08)

write.csv(coloc_output_08, 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_STARNET_merged_eQTLs_HPP4_08_P_5e_8_Renamed_pheno_renamed_tissue.csv', row.names = FALSE)
