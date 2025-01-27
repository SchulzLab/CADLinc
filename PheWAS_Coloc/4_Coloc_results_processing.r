library(dplyr)
library(ggplot2)
library(plyr)
library(tidyr)
library(tidyverse)
library(extrafont)
library(scales)

# I. Get GTEx effect eQTL alleles (didn't include that info in the beginning of the analysis) 
coloc_GTEx <- read.csv("results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_eQTLs_HPP4_06_P_5e_8.csv")

snps_to_extract <- unique(coloc_GTEx[, c("SNP", "chr_pos_hg38")])

# Get the list of relevant GTEx eQTL files
# Select only relevant files and tissues (13) 
# Directory path
directory <- "GTEx_Analysis_v8_eQTL/"

# List all files in the directory
files <- list.files(directory, full.names = TRUE)

# Filter only files with "signif_variant_gene_pairs" in their names
all_variant_gene_pairs_files <- files[grep("signif_variant_gene_pairs", files)]

# Further filter the files to include only those with the desired substrings
desired_substrings <- c('Artery_Aorta', 'Whole_Blood', 'Artery_Coronary', 'Liver', 'Adipose_Subcutaneous', 'Muscle_Skeletal', 'Artery_Tibial', 'Adipose_Visceral_Omentum', 'Heart_Atrial_Appendage',
                  'Heart_Left_Ventricle', 'Nerve_Tibial', 'Kidney_Cortex', 'Spleen')

pattern <- paste(desired_substrings, collapse = "|")
variant_gene_pairs_files <- all_variant_gene_pairs_files[grep(pattern, all_variant_gene_pairs_files)]

tissue_names <- gsub(".*/([^/]+)\\.v8.*", "\\1", variant_gene_pairs_files)

# Run a cycle to find eQTLs from GTEx data for all SNPs, all tissues
# Initialize an empty dataframe to store the final merged results
results <- data.frame()

# Iterate over dataset IDs from datasets_filtered dataframe
for (file in variant_gene_pairs_files) {
  
  # Read the file
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Add the column with only chr and bp position
  data$chr_pos_hg38 <- sub("^(.*)_[^_]*$", "\\1", data$variant_id)
  data$chr_pos_hg38 <- sub("^(.*)_[^_]*$", "\\1", data$chr_pos_hg38)
  data$chr_pos_hg38 <- sub("^(.*)_[^_]*$", "\\1", data$chr_pos_hg38)
  data$chr_pos_hg38 <- gsub("chr", "", data$chr_pos_hg38) 
  data$chr_pos_hg38 <- gsub("_", ":", data$chr_pos_hg38)
  
  data$gene_ensembl <- sub("\\.\\d+", "", data$gene_id)
  
  # Get the tissue name
  tissue_name <- gsub(".*/([^/]+)\\.v8.*", "\\1", file)
  # Add the column with the tissue
  data$tissue <- tissue_name
  
  # Find the overlap in dataframes based on SNP hg38 position
  gene_eQTL_merged <- merge(snps_to_extract, data, by = "chr_pos_hg38")
  
  results <- rbind(results, gene_eQTL_merged)
  
  # Print a message indicating that the tissue has been processed
  cat("Tissue", tissue_name, "has been processed.\n")

}

results_subset <- (results[, c("SNP", "chr_pos_hg38", 'gene_id', 'slope', 'pval_nominal', 'tissue', 'variant_id')])

merged_df <- coloc_GTEx %>%
  inner_join(results_subset, by = c("SNP", "chr_pos_hg38", "gene_id", "slope", "pval_nominal", "tissue"))

# Separate the variant_id column into multiple columns
merged_df <- merged_df %>%
  separate(variant_id, into = c("chr_eQTL", "pos_eQTL", "A2_eQTL", "A_effect_eQTL", "b38"), sep = "_", remove = FALSE)

merged_df$b38 <- NULL
merged_df$chr_eQTL <- NULL
merged_df$pos_eQTL <- NULL

write.csv(merged_df, "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_eQTLs_HPP4_06_P_5e_8_with_eQTL_alleles.csv", row.names = FALSE)

# II. Process coloc output files
# Data preparation
coloc_output_GTEx <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_eQTLs_HPP4_06_P_5e_8_with_eQTL_alleles.csv')
coloc_output_GTEx$Study <- 'GTEx'

coloc_output_STARNET <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_STARNET_eQTLs_HPP4_06_P_5e_8.csv')
coloc_output_STARNET$Study <- 'STARNET'

STARNET_eQTL_alleles <- read.table('genotype_case_snp_imputation_info.txt', sep = '\t', header = TRUE)
STARNET_eQTL_alleles <- STARNET_eQTL_alleles[, c("marker_id", "a1", 'a2')]
names(STARNET_eQTL_alleles)[names(STARNET_eQTL_alleles) == "marker_id"] <- "SNP"
names(STARNET_eQTL_alleles)[names(STARNET_eQTL_alleles) == "a1"] <- "A_effect_eQTL"
names(STARNET_eQTL_alleles)[names(STARNET_eQTL_alleles) == "a2"] <- "A2_eQTL"

coloc_output_STARNET <- merge(coloc_output_STARNET, STARNET_eQTL_alleles, by = c("SNP"))

# Merge 2 coloc outputs and select only similar columns
# Make a subset
coloc_output_GTEx <- coloc_output_GTEx[, c("SNP", "position", "SNP.PP.H4", "gene_ensembl", "chr", "position_hg19", "chr_bp_position_hg19", "slope", "pval_nominal", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", "SE", "sample_size", "tissue", "Trait", "test", "Gene.name",  "Biotype", "Sample_size_eQTL", "Study", "A_effect_eQTL", "A2_eQTL")]

coloc_output_STARNET <- coloc_output_STARNET[, c("SNP", "position", "SNP.PP.H4", "gene_ensembl", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "p_value_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", "SE", "sample_size", "tissue", "Trait", "test", "Gene.name",  "Biotype", "Sample_size_eQTL", "Study", "A_effect_eQTL", "A2_eQTL")]

# Rename some columns
names(coloc_output_GTEx)[names(coloc_output_GTEx) == "slope"] <- "beta_eQTL"
names(coloc_output_GTEx)[names(coloc_output_GTEx) == "pval_nominal"] <- "p_value_eQTL"

coloc_output <- rbind(coloc_output_GTEx, coloc_output_STARNET)

# Add new column with 2 types in biotype
coloc_output <- coloc_output %>%
  mutate(Biotype_binary = ifelse(Biotype == "protein_coding", "protein-coding", "non-coding"))
table(coloc_output$Biotype_binary)

# Rename traits for plotting
coloc_output$Trait_labels <- coloc_output$Trait

coloc_output$Trait_labels <- revalue(coloc_output$Trait_labels,
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

head(coloc_output)

# Rename Tissue for plotting, merge them in froups
coloc_output$tissue_labels <- coloc_output$tissue

coloc_output$tissue_labels <- revalue(coloc_output$tissue_labels,
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
                                     "Nerve_Tibial" = "Tibial nerve"
                                     ))

head(coloc_output)

write.csv(coloc_output, 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_STARNET_merged_eQTLs_HPP4_06_P_5e_8_with_eQTL_alleles_Renamed_pheno_renamed_tissue.csv', row.names = FALSE)

