library(tidyverse)   
library(data.table)   
library(stringr)     
library(readr)        
library(plyr)         
library(rtracklayer)

# File with 1580 candidate CAD genes
JointGene_Bodies <- read.table('data/230524_FixDC_GeneTable.txt', sep = '\t', header = TRUE)
head(JointGene_Bodies)

# I. STARNET
# Get eQTL SNPs

# Get the list of all relevant STARNET eQTL files
# Directory path
directory <- "STARNET_eQTL_FDR01/"

# List all files in the directory
files <- list.files(directory, full.names = TRUE)

# Filter only files with " FDR01.tbl" in their names
starnet_files <-list.files(directory, pattern = "FDR01\\.tbl$", full.names = TRUE)
tissue_names <- sub(".*/([^_/]+).*", "\\1", starnet_files)

# Run a cycle to find eQTLs from GTEx data for all SNPs, all tissues
# Initialize an empty dataframe to store the final merged results
results <- data.frame()

# Iterate over dataset IDs from datasets_filtered dataframe
for (file in starnet_files) {
  
  # Read the file
  data <- read.table(file, header = TRUE, sep = " ")
  
  # Get the tissue name
  tissue_name <- sub(".*/([^_/]+).*", "\\1", file)
  # Add the column with the tissue
  data$tissue <- tissue_name
  data[, tissue_name] <- NULL
  
  # Find the overlap in dataframes based on the gene id
  gene_eQTL_merged <- merge(JointGene_Bodies, data, by.x = "gene_ensembl", by.y = "gene")
  
  results <- rbind(results, gene_eQTL_merged)
  
  # Print a message indicating that the tissue has been processed
  cat("Tissue", tissue_name, "has been processed.\n")

}

write.csv(results, 'results/gene_eQTL_GWAS_Joint_GeneBase/All_eQTL_STARNET_JointGenes_SNEEP_GATESGeneBodies.csv', row.names = FALSE)


# II. GTEx 
# 1. 
# Get eQTL SNPs

# Get the list of all relevant GTEx eQTL files
# Directory path
directory <- "GTEx_Analysis_v8_eQTL/"

# List all files in the directory
files <- list.files(directory, full.names = TRUE)

# Filter only files with "signif_variant_gene_pairs" in their names
variant_gene_pairs_files <- files[grep("signif_variant_gene_pairs", files)]
tissue_names <- gsub(".*/([^/]+)\\.v8.*", "\\1", variant_gene_pairs_files)

# Run a cycle to find eQTLs from GTEx data for all SNPs, all tissues
# Initialize an empty dataframe to store the final merged results
results <- data.frame()

# Iterate over dataset IDs from datasets_filtered dataframe
for (file in variant_gene_pairs_files) {
  
  # Read the file
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Add the column with only chr and bp position
  data$chr_bp <- sub("^(.*)_[^_]*$", "\\1", data$variant_id)
  data$chr_bp <- sub("^(.*)_[^_]*$", "\\1", data$chr_bp)
  data$chr_bp <- sub("^(.*)_[^_]*$", "\\1", data$chr_bp)
  
  data$gene_ensembl <- sub("\\.\\d+", "", data$gene_id)
  
  # Get the tissue name
  tissue_name <- gsub(".*/([^/]+)\\.v8.*", "\\1", file)
  # Add the column with the tissue
  data$tissue <- tissue_name
  
  # Find the overlap in dataframes based on the gene id
  gene_eQTL_merged <- merge(JointGene_Bodies, data, by = "gene_ensembl")
  
  results <- rbind(results, gene_eQTL_merged)
  
  # Print a message indicating that the tissue has been processed
  cat("Tissue", tissue_name, "has been processed.\n")

}

results$chr_bp <- NULL
write.csv(results, 'results/gene_eQTL_GWAS_Joint_GeneBase/All_eQTL_GTEx_JointGenes_SNEEP_GATESGeneBodies.csv', row.names = FALSE)


# 2.
# Process and prepare eQTL SNPs data (add hg19 position and rsID) before PheWAS

# 2.1. Get hg19 position

# Get the "chr" and "position_hg38" columns for GTEx data
GTEx_eQTL$chr_bp <- sub("^(.*)_[^_]*$", "\\1", GTEx_eQTL$variant_id)
GTEx_eQTL$chr_bp <- sub("^(.*)_[^_]*$", "\\1", GTEx_eQTL$chr_bp)
GTEx_eQTL$chr_bp <- sub("^(.*)_[^_]*$", "\\1", GTEx_eQTL$chr_bp)

# Get hg38 position
GTEx_eQTL$position_hg38 <- sapply(strsplit(GTEx_eQTL$chr_bp, "_"), function(x) x[[2]])

# Get chr
GTEx_eQTL$chr <- sub("_.*", "", GTEx_eQTL$chr_bp)
# Extract numeric part
GTEx_eQTL$chr <- as.numeric(sub("chr", "", GTEx_eQTL$chr))

unique_SNPs_positions <- unique(GTEx_eQTL[, c("chr", "position_hg38")])
write.table(unique_SNPs_positions[, c("chr", "position_hg38")], file = 'results/gene_eQTL_GWAS_Joint_GeneBase/pars_rsID_SNPs_Unique_eQTL_GTEx_JointGenes_SNEEP_GATESGeneBodies.txt', sep = "\t", row.names = FALSE, col.names = FALSE)

# Prepare bed files to extract hg19 position with liftover tool (online uploading the file)
# Get a new column with chr like chr1
unique_SNPs_positions$chr_hg38 <- paste0('chr', unique_SNPs_positions$chr)

# Create a character vector with BED format
bed_lines <- paste0(unique_SNPs_positions$chr_hg38, "\t", unique_SNPs_positions$position_hg38, "\t", unique_SNPs_positions$position_hg38)

# Write the BED lines to a file
writeLines(bed_lines, "results/gene_eQTL_GWAS_Joint_GeneBase/hg38_All_eQTL_GTEx_JointGenes_SNEEP_GATESGeneBodies.bed")

# Go to LiftOver online tool
# https://genome.ucsc.edu/cgi-bin/hgLiftOver

# Read the table with hg19 position
snps_hg19 <- read.table('results/gene_eQTL_GWAS_Joint_GeneBase/hg19_All_eQTL_GTEx_JointGenes_SNEEP_GATESGeneBodies.bed')
head(snps_hg19)

# Add the column for chr_bp hg38
GTEx_eQTL$chr_pos_hg38 <- paste(GTEx_eQTL$chr, GTEx_eQTL$position_hg38, sep = ":")
head(GTEx_eQTL)

# Remove 'chr' part in the chr number
snps_hg19$V1 <- sub("chr", "", snps_hg19$V1)
colnames(snps_hg19)[colnames(snps_hg19) == "V1"] <- "chr_hg19"
colnames(snps_hg19)[colnames(snps_hg19) == "V2"] <- "position_hg19"

# Add the column chr_bp hg19 to the df from liftover
# Remove the part after ":" and before "-"
snps_hg19$chr_pos_hg38 <- gsub(":.*?-(\\d+)", ":\\1", snps_hg19$V4)
# Remove 'chr' part from the values in chr_pos_hg38
snps_hg19$chr_pos_hg38 <- sub('^chr', '', snps_hg19$chr_pos_hg38)
head(snps_hg19)

# Add the column with hg19 from snps_hg19 to the GTEx_eQTL
merged_hg19_GTEx_eQTL <- merge(snps_hg19, GTEx_eQTL, by = "chr_pos_hg38", all = TRUE)
merged_hg19_GTEx_eQTL <- unique(merged_hg19_GTEx_eQTL)
merged_hg19_GTEx_eQTL$V3 <-  NULL
merged_hg19_GTEx_eQTL$V4 <-  NULL
merged_hg19_GTEx_eQTL$V5 <-  NULL

# Add the column chr:position for hg19
merged_hg19_GTEx_eQTL$chr_bp_position_hg19 <- paste(merged_hg19_GTEx_eQTL$chr_hg19, merged_hg19_GTEx_eQTL$position_hg19, sep = ":")

write.csv(merged_hg19_GTEx_eQTL, 'results/gene_eQTL_GWAS_Joint_GeneBase/hg19_hg38_All_eQTL_GTEx_Genes_JointGenes_SNEEP_GATESGeneBodies.csv', row.names = FALSE)

# 2.2. Get rsIDs

# Prepare the file for rsID parsing from dbSNP hg38 file
snps_pars <- read.table('results/gene_eQTL_GWAS_Joint_GeneBase/pars_rsID_SNPs_Unique_eQTL_GTEx_JointGenes_SNEEP_GATESGeneBodies.txt', header = FALSE, stringsAsFactors = FALSE, fill = TRUE, row.names = NULL)

snps_pars_prep <- data.frame(
  V1 = snps_pars$V1,
  Second = snps_pars$V2 - 1,
  Third = snps_pars$V2
)

# Sort rows
snps_pars_prep <- snps_pars_prep[order(snps_pars_prep$V1, snps_pars_prep$Second), ]

# Print the first few rows of the new dataframe
head(snps_pars_prep)

# Define the order of chromosome names in V1
chromosome_order <- c('1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9')

# Sort the rows
snps_pars_prep <- snps_pars_prep[order(factor(snps_pars_prep$V1, levels = chromosome_order), snps_pars_prep$Second), ]

write.table(snps_pars_prep, file = "results/gene_eQTL_GWAS_Joint_GeneBase/snps_pars_genes_natur_sort.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Run the command
# bedmap --echo --echo-map-id-uniq snps_pars_genes_natur_sort.bed sorted_common_all_20180418.starch > parsed_430k_snps_with_rsID_natur_sort.bed

# Process the file
file_path <- "results/gene_eQTL_GWAS_Joint_GeneBase/parsed_430k_snps_with_rsID_natur_sort.bed"

# Read the file line by line
lines <- readLines(file_path)

# Split each line by tabs
data_elements <- strsplit(lines, "\t")

# Extract chromosome, start, and combined end|rsID columns
chromosome <- as.numeric(sapply(data_elements, `[[`, 1))
start <- as.numeric(sapply(data_elements, `[[`, 2))

# Split the third element (end|rsID) by "|"
end_rsID <- sapply(data_elements, `[[`, 3)
end_rsID_split <- strsplit(end_rsID, "\\|")

# Extract the end position and rsID from each split
end <- as.numeric(sapply(end_rsID_split, `[[`, 1))
rsID <- sapply(end_rsID_split, function(x) if(length(x) > 1) x[2] else NA)

# Combine the extracted columns into a data frame
all_parsed_snps_rsIDs <- data.frame(chromosome, start, end, rsID)

# Print the first few rows to check if it was parsed correctly
head(all_parsed_snps_rsIDs)

# Split rows with multiple rsIDs into separate rows
df_split <- all_parsed_snps_rsIDs %>%
  separate_rows(rsID, sep = ";") %>%
  mutate(rsID = ifelse(rsID == "", NA, rsID))  # Replace empty strings with NA

# Fill NA values in rsID column with the previous non-NA value in the same group
df_filled <- df_split %>%
  group_by(chromosome, start, end) %>%
  fill(rsID)

# Remove duplicate rows
all_parsed_snps_rsIDs_final <- df_filled %>%
  distinct()

head(all_parsed_snps_rsIDs_final)  

# Add the column with hg38 position
all_parsed_snps_rsIDs_final$chr_pos_hg38 <- paste(all_parsed_snps_rsIDs_final$chromosome, all_parsed_snps_rsIDs_final$end, sep = ":")

write.csv(all_parsed_snps_rsIDs_final, 'results/gene_eQTL_GWAS_Joint_GeneBase/all_parsed_430k_snps_rsID_final_processed_hg38.csv', row.names = FALSE)

# Add rsIDs to the the GTEx eQTL file
GTEX_hg19_38_rsID <- merge(merged_hg19_GTEx_eQTL, all_parsed_snps_rsIDs_final, by = "chr_pos_hg38", all = TRUE)

GTEX_hg19_38_rsID <- GTEX_hg19_38_rsID[complete.cases(GTEX_hg19_38_rsID$rsID), ]
names(GTEX_hg19_38_rsID)[names(GTEX_hg19_38_rsID) == "rsID"] <- "SNP"

write.csv(GTEX_hg19_38_rsID, 'results/gene_eQTL_GWAS_Joint_GeneBase/hg19_hg38_rsID_All_eQTL_GTEx_Genes_JointGenes_SNEEP_GATESGeneBodies.csv', row.names = FALSE)
