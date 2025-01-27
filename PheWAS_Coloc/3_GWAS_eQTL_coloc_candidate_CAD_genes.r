library(coloc)
library(doParallel)
library(foreach)
library(dplyr)
library(plyr)

# I. GTEx
# 1. Prepare the files for coloc
GTEx_merged <- read.csv("results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Final_merged_tables/PheWAS_results_GTEx_eQTLs_FINAL_Genes_Joint_GeneBase_5e_5.csv")
head(GTEx_merged)

# Select relevant tissue types
GTEx_merged <- GTEx_merged[GTEx_merged$tissue %in% c('Artery_Aorta', 'Whole_Blood', 'Artery_Coronary', 'Liver', 'Adipose_Subcutaneous', 'Muscle_Skeletal', 'Artery_Tibial', 'Adipose_Visceral_Omentum', 'Heart_Atrial_Appendage',
                  'Heart_Left_Ventricle', 'Nerve_tibial', 'Kidney_Cortex', 'Spleen'), ]
GTEx_merged <- GTEx_merged[GTEx_merged$p_value_GWAS < 5e-8, ]

GTEx_merged$position_hg19 <- sub(".*:", "", GTEx_merged$chr_bp_position_hg19)

write.csv(GTEx_merged, 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/PheWAS_GTEx_eQTLs_13_tissues_FINAL_5e_8.csv', row.names = FALSE)

# For 15 SNPs the hg19 position is missing
rows_with_NA_in_hp19 <- GTEx_merged[is.na(GTEx_merged$position_hg19), ]
table(rows_with_NA_in_hp19$Trait)

# For those SNPs manually add the position
GTEx_merged[GTEx_merged$SNP == "rs1137571" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("19:45219011", 45219011)

GTEx_merged[GTEx_merged$SNP == "rs74862042" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("19:45217784", 45217784)

GTEx_merged[GTEx_merged$SNP == "rs1810741" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("20:34109715", 34109715)

GTEx_merged[GTEx_merged$SNP == "rs224437" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("20:34154371", 34154371)

GTEx_merged[GTEx_merged$SNP == "rs28680494" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("22:42528851", 42528851)

GTEx_merged[GTEx_merged$SNP == "rs35132383" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("1:150124780", 150124780)

GTEx_merged[GTEx_merged$SNP == "rs1815302" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("1:150257933", 150257933)

GTEx_merged[GTEx_merged$SNP == "rs9849509" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("3:48483432", 48483432)

GTEx_merged[GTEx_merged$SNP == "rs9876781" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("3:48487338", 48487338)

GTEx_merged[GTEx_merged$SNP == "rs551154" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("9:136197278", 136197278)

GTEx_merged[GTEx_merged$SNP == "rs642059" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("9:136195734", 136195734)

GTEx_merged[GTEx_merged$SNP == "rs549443" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("9:136135237", 136135237)

GTEx_merged[GTEx_merged$SNP == "rs5011221" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("1:149772497", 149772497)

GTEx_merged[GTEx_merged$SNP == "rs576123" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("9:136144308", 136144308)

GTEx_merged[GTEx_merged$SNP == "rs61285422" & is.na(GTEx_merged$position_hg19), c("chr_bp_position_hg19", "position_hg19")] <- list("19:45219903", 45219903)

# 2. Run coloc
path_for_MAF <- 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Intermediate_files_for_MAF_extraction'

# Make a function
# Final cycle: all traits and all genes one by one by tissue
# Define a function to process each trait
process_trait <- function(trait_name) {
    library(dplyr)
    library(plyr)
    library(coloc)
  
    # Make subset based on trait value
    GWAS_eQTL_trait <- GTEx_merged[GTEx_merged$Trait == trait_name, ]
    
    # 1) 
    # Add MAF info using PLINK
    # Step 1: Write SNP IDs to a text file
    write.table(GWAS_eQTL_trait$SNP, paste0(path_for_MAF, '/', trait_name, "/snps_to_extract.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Step 2: Execute PLINK command to calculate allele frequencies
    plink_command <- paste("/PLINK_tool/plink",
                       "--bfile /data/PLINK_reference/Merged_files/merged_bfile",
                       "--extract", paste0(path_for_MAF, '/', trait_name, "/snps_to_extract.txt"),
                       "--freq --out", paste0(path_for_MAF, '/', trait_name, "/plink_MAF_output"))
    
    system(plink_command)
    
    # Step 3: Parse PLINK output to extract MAF values
    plink_output <- read.table(paste0(path_for_MAF, '/', trait_name, "/plink_MAF_output.frq"), header = TRUE)
    MAF_values <- plink_output$MAF
    
    # Step 4: Merge MAF values back into the original data frame
    GWAS_eQTL_trait$MAF <- MAF_values[match(GWAS_eQTL_trait$SNP, plink_output$SNP)]
    
    # Select only SNPs with available MAF values (present in EUR population)
    GWAS_eQTL_trait <- GWAS_eQTL_trait[complete.cases(GWAS_eQTL_trait$MAF), ]
    
    # 2)
    # Add the sample size for eQTL part
    GWAS_eQTL_trait <- GWAS_eQTL_trait %>%
  mutate(
    Sample_size_eQTL = case_when(
      tissue == "Adipose_Subcutaneous" ~ 581,
      tissue == "Adipose_Visceral_Omentum" ~ 469,
      tissue == "Artery_Aorta" ~ 387,
      tissue == "Artery_Coronary" ~ 213,
      tissue == "Artery_Tibial" ~ 584,
      tissue == "Heart_Atrial_Appendage" ~ 372,
      tissue == "Heart_Left_Ventricle" ~ 386,
      tissue == "Kidney_Cortex" ~ 73,
      tissue == "Liver" ~ 208,
      tissue == "Muscle_Skeletal" ~ 706,
      tissue == "Nerve_Tibial" ~ 532,
      tissue == "Spleen" ~ 227,
      tissue == "Whole_Blood" ~ 670,
    )
  )
     
    # Create an empty dataframe to store results for the current trait
    trait_results <- data.frame()
  
    # Iterate over genes in the current trait
    unique_genes <- unique(GWAS_eQTL_trait$gene_ensembl)
    for (gene in unique_genes) {
      
      # Filter data for the current gene
      one_gene <- GWAS_eQTL_trait[GWAS_eQTL_trait$gene_ensembl == gene, ]
      
        # Iterate over unique tissue values
      unique_tissues <- unique(one_gene$tissue)
      for (tissue in unique_tissues) {
        # Filter data for the current tissue
        one_tissue <- one_gene[one_gene$tissue == tissue, ]
    
    # !!! There are some duplicated SNPs in eQTL data (the same SNP is mentioned a couple of times for the same gene in the same tissue but has different values in the variables)
    # Leave the row with the lowest eQTL p-value for each SNP
    # Convert 'p_value_eQTL' to numeric if it's not already
    
    names(one_tissue)[names(one_tissue) == "pval_nominal"] <- "p_value_eQTL"    
    
    one_tissue$p_value_eQTL <- as.numeric(one_tissue$p_value_eQTL)
    one_tissue <- one_tissue %>%
      group_by(SNP) %>%
      filter(p_value_eQTL == min(p_value_eQTL)) %>%
      ungroup()
    
      # Filter then by GWAS p-value (if the duplicates are from eQTL part)
      one_tissue <- one_tissue %>%
      group_by(SNP) %>%
      filter(p_value_GWAS == min(p_value_GWAS)) %>%
      ungroup()
    
       # Check for completely duplicated rows and remove them
      one_tissue <- one_tissue %>%
      distinct()  # Remove completely duplicated rows
      
      # Prepare 2 files: with GWAS and eQTL data 
      # 1.
      # Make a subset of eQTL part for the tool
      eQTL_trait <- one_tissue[, c("SNP", 'position_hg19', "chr", 'slope', 'p_value_eQTL', 'gene_ensembl', 'Gene.name', 'test', 'MAF', 'Sample_size_eQTL')]
      
      eQTL_trait$type <- 'quant'

      df1=eQTL_trait
      colnames(df1) <- c("snp", "position", "chr", "beta", 'pvalues', 'ensembl', 'Gene.name', 'test', 'MAF', 'N', "type")
      df1 <- as.list(df1)

      
      # 2.
      # Make a subset of GWAS part for the tool
      GWAS_trait <- one_tissue[, c("SNP", 'position_hg19', "chr", 'beta_GWAS', 'p_value_GWAS', 'gene_ensembl', 'Gene.name', 'test', 'MAF','sample_size', 'SE')]
      
      
      GWAS_trait <- GWAS_trait %>%
  mutate(
    type = case_when(
      trait_name %in% c("NAFLD", "pad_eversmoker", "pad_nodiabetes", 'Atrial_fibrillation', 'Heart_Failure', 'pad_diabetes', 'pad_neversmoker', 'pad_primary', 'Type_2_diabetes') ~ "cc",
      trait_name %in% c("Basophil_count", "Eosinophil_count", "Fasting_insulin", 'HDL_cholesterol', 'Lymphocyte_count', 'Non-HDL_cholesterol', 'Platelet_count', 'Red_blood_cell_count', 'Systolic_blood_pressure', 'Triglycerides', 'Waist-hip_ratio', 'Waist-hip_ratio_adj_BMI', 'Body_mass_index', 'Diastolic_blood_pressure', 'Fasting_glucose', 'HbA1c', 'LDL_cholesterol', 'Monocyte_count', 'Neutrophil_count', 'Pulse_pressure', 'Resting_heart_rate', 'Total_cholesterol', 'Two_hour_glucose', 'White_blood_cell_count') ~ "quant",
    )
  )
      
      GWAS_trait$varbeta <- GWAS_trait$SE^2
      
      df2=GWAS_trait
      colnames(df2) <- c("snp", "position", "chr", "beta",'pvalues', 'ensembl', 'Gene.name', 'test', 'MAF', "N", 'SE', "type", "varbeta")
      df2 <- as.list(df2)
      
      
      # Convert data frames to lists
      df1_list <- list(
        snp = df1$snp,
        position = df1$position,
        chr = df1$chr,
        beta = df1$beta,
        pvalues = df1$pvalues,
        ensembl = df1$ensembl,
        Gene.name = df1$Gene.name,
        test = df1$test,
        MAF = df1$MAF,
        N = df1$N,
        type = unique(df1$type)[1]  # Ensure type is a single value
      )
      
      df2_list <- list(
        snp = df2$snp,
        position = df2$position,
        chr = df2$chr,
        beta = df2$beta,
        pvalues = df2$pvalues,
        ensembl = df2$ensembl,
        Gene.name = df2$Gene.name,
        test = df2$test,
        MAF = df2$MAF,
        N = df2$N,
        SE = df2$SE,
        type = unique(df2$type)[1],  # Ensure type is a single value
        varbeta = df2$varbeta
      )
      
      # Run coloc
      
      my.res <- coloc.abf(dataset1 = df1_list, dataset2 = df2_list)
      
      # Get results
      results_subset <- subset(my.res$results)
      colnames(results_subset)[colnames(results_subset) == "snp"] <- "SNP"
      selected_snps <- results_subset$SNP
      selected_GWAS_eQTL <- GWAS_eQTL_trait[GWAS_eQTL_trait$SNP %in% selected_snps, ]
      
      # Merge results
      coloc_results <- merge(results_subset, selected_GWAS_eQTL, by = "SNP", all = TRUE)
      
      # Append to the trait_results data frame
      trait_results <- rbind(trait_results, coloc_results)
      
    }
      
    cat("Processed:", gene, ", ", trait_name, ", Tissue:", tissue, "\n")
    
  }
    
  # Save results for the current trait to a CSV file
  write.csv(trait_results, file = paste0('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/GTEx_results/GTEx_coloc_results_', trait_name, '.csv'), row.names = FALSE)
    
  cat("Processed trait:", trait_name, "\n")
    
  # Make sure to return the trait_results dataframe
  return(trait_results)
}

# Select phenotypes
unique_traits <- c("Body_mass_index", "Waist-hip_ratio_adj_BMI", "Waist-hip_ratio",
                   "LDL_cholesterol", "HDL_cholesterol", "Total_cholesterol", "Non-HDL_cholesterol", "Triglycerides",
                   "Fasting_insulin", "Two_hour_glucose", "Fasting_glucose", "HbA1c",
                   "Diastolic_blood_pressure", "Pulse_pressure", "Systolic_blood_pressure", 
                   "Neutrophil_count", "White_blood_cell_count", "Lymphocyte_count", "Basophil_count", "Monocyte_count", "Platelet_count", "Red_blood_cell_count", "Eosinophil_count",
                  "Type_2_diabetes", "Atrial_fibrillation", "Heart_Failure", "pad_eversmoker", "pad_nodiabetes", "pad_primary",  "NAFLD", 'Resting_heart_rate', "CAD")

# Parallel calculation
as.integer(system("nproc --all", intern = TRUE))

# Set up parallel backend with desired number of cores
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Iterate over unique trait values in parallel
trait_results_list <- foreach(trait_name = unique_traits, .combine = rbind) %dopar% {
    process_trait(trait_name)
}

# Close the parallel backend
stopCluster(cl)

# 3. Merge the coloc results
# make a subset based on HPP4 >= 0.6
# SNP.PP.H4 column
# Set your folder path
folder_path <- "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/GTEx_results/"

# Get a list of CSV files in the folder
# GTEx
csv_files <- list.files(path = folder_path, pattern = '^GTEx', full.names = TRUE)

# Initialize an empty data frame to store the merged data
merged_data <- data.frame() 
 
# Loop through the CSV files, read them, and row-bind to the merged_data
for (csv_file in csv_files) {
  data <- read.csv(csv_file, header = TRUE)
  data <- subset(data, SNP.PP.H4 >= 0.6)
  merged_data <- rbind(merged_data, data)
  
  # Print the file being processed
  cat("Processed file:", csv_file, "\n")
}  

write.csv(merged_data, file = "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_eQTLs_HPP4_06_P_5e_8.csv", row.names = FALSE)


# II. STARNET data
# 1. Prepare the files for coloc
STARNET_merged <- read.csv("results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Final_merged_tables/PheWAS_results_STARNET_eQTLs_FINAL_Genes_Joint_GeneBase_5e_5.csv")
head(STARNET_merged)

# All tissue types are already relevant
STARNET_merged <- STARNET_merged[STARNET_merged$p_value_GWAS < 5e-8, ]
write.csv(STARNET_merged, 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/PheWAS_STARNET_eQTLs_5e_8.csv', row.names = FALSE)  

# 2. Run coloc
path_for_MAF <- 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Intermediate_files_for_MAF_extraction/STARNET' 

# Make a function
# Define a function to process each trait
process_trait <- function(trait_name) {
    library(dplyr)
    library(plyr)
    library(coloc)
  
    # Make subset based on trait value
    GWAS_eQTL_trait <- STARNET_merged[STARNET_merged$Trait == trait_name, ]
    
    # 1) 
    # Add MAF info using PLINK
    # Step 1: Write SNP IDs to a text file
    write.table(GWAS_eQTL_trait$SNP, paste0(path_for_MAF, '/', trait_name, "/snps_to_extract.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Step 2: Execute PLINK command to calculate allele frequencies
    plink_command <- paste("/PLINK_tool/plink",
                       "--bfile /data/PLINK_reference/Merged_files/merged_bfile",
                       "--extract", paste0(path_for_MAF, '/', trait_name, "/snps_to_extract.txt"),
                       "--freq --out", paste0(path_for_MAF, '/', trait_name, "/plink_MAF_output"))
    
    system(plink_command)
    
    # Step 3: Parse PLINK output to extract MAF values
    plink_output <- read.table(paste0(path_for_MAF, '/', trait_name, "/plink_MAF_output.frq"), header = TRUE)
    MAF_values <- plink_output$MAF
    
    # Step 4: Merge MAF values back into the original data frame
    GWAS_eQTL_trait$MAF <- MAF_values[match(GWAS_eQTL_trait$SNP, plink_output$SNP)]
    
    # Select only SNPs with available MAF values (present in EUR population)
    GWAS_eQTL_trait <- GWAS_eQTL_trait[complete.cases(GWAS_eQTL_trait$MAF), ]
    
    # 2)
    # Add the sample size for eQTL part
    GWAS_eQTL_trait <- GWAS_eQTL_trait %>%
  mutate(
    Sample_size_eQTL = case_when(
      tissue == "AOR" ~ 600,
      tissue == "BLD" ~ 600,
      tissue == "LIV" ~ 600,
      tissue == "MAM" ~ 600,
      tissue == "SF" ~ 600,
      tissue == "SKLM" ~ 600,
      tissue == "VAF" ~ 600
    )
  )
     
    # Create an empty dataframe to store results for the current trait
    trait_results <- data.frame()
  
    # Iterate over genes in the current trait
    unique_genes <- unique(GWAS_eQTL_trait$gene_ensembl)
    for (gene in unique_genes) {
      
      # Filter data for the current gene
      one_gene <- GWAS_eQTL_trait[GWAS_eQTL_trait$gene_ensembl == gene, ]
      
        # Iterate over unique tissue values
      unique_tissues <- unique(one_gene$tissue)
      for (tissue in unique_tissues) {
        # Filter data for the current tissue
        one_tissue <- one_gene[one_gene$tissue == tissue, ]
    
    # !!! There are some duplicated SNPs in eQTL data (the same SNP is mentioned a couple of times for the same gene in the same tissue but has different values in the variables)
    # Leave the row with the lowest eQTL p-value for each SNP
    # Convert 'p_value_eQTL' to numeric if it's not already
    
    names(one_tissue)[names(one_tissue) == "pval_nominal"] <- "p_value_eQTL"    
    
    one_tissue$p_value_eQTL <- as.numeric(one_tissue$p_value_eQTL)
    one_tissue <- one_tissue %>%
      group_by(SNP) %>%
      filter(p_value_eQTL == min(p_value_eQTL)) %>%
      ungroup()
    
      # Filter then by GWAS p-value (if the duplicates are from eQTL part)
      one_tissue <- one_tissue %>%
      group_by(SNP) %>%
      filter(p_value_GWAS == min(p_value_GWAS)) %>%
      ungroup()
    
       # Check for completely duplicated rows and remove them
      one_tissue <- one_tissue %>%
      distinct()  # Remove completely duplicated rows
      
      # Prepare 2 files: with GWAS and eQTL data 
      # 1.
      # Make a subset of eQTL part for the tool
      eQTL_trait <- one_tissue[, c("SNP", 'position_hg19', "chr", 'beta_eQTL', 'p_value_eQTL', 'gene_ensembl', 'Gene.name', 'test', 'MAF', 'Sample_size_eQTL')]
      
      eQTL_trait$type <- 'quant'

      df1=eQTL_trait
      colnames(df1) <- c("snp", "position", "chr", "beta", 'pvalues', 'ensembl', 'Gene.name', 'test', 'MAF', 'N', "type")
      df1 <- as.list(df1)

      
      # 2.
      # Make a subset of GWAS part for the tool
      GWAS_trait <- one_tissue[, c("SNP", 'position_hg19', "chr", 'beta_GWAS', 'p_value_GWAS', 'gene_ensembl', 'Gene.name', 'test', 'MAF','sample_size', 'SE')]
      
      
      GWAS_trait <- GWAS_trait %>%
  mutate(
    type = case_when(
      trait_name %in% c("NAFLD", "pad_eversmoker", "pad_nodiabetes", 'Atrial_fibrillation', 'Heart_Failure', 'pad_diabetes', 'pad_neversmoker', 'pad_primary', 'Type_2_diabetes') ~ "cc",
      trait_name %in% c("Basophil_count", "Eosinophil_count", "Fasting_insulin", 'HDL_cholesterol', 'Lymphocyte_count', 'Non-HDL_cholesterol', 'Platelet_count', 'Red_blood_cell_count', 'Systolic_blood_pressure', 'Triglycerides', 'Waist-hip_ratio', 'Waist-hip_ratio_adj_BMI', 'Body_mass_index', 'Diastolic_blood_pressure', 'Fasting_glucose', 'HbA1c', 'LDL_cholesterol', 'Monocyte_count', 'Neutrophil_count', 'Pulse_pressure', 'Resting_heart_rate', 'Total_cholesterol', 'Two_hour_glucose', 'White_blood_cell_count') ~ "quant",
    )
  )
      
      GWAS_trait$varbeta <- GWAS_trait$SE^2
      
      df2=GWAS_trait
      colnames(df2) <- c("snp", "position", "chr", "beta",'pvalues', 'ensembl', 'Gene.name', 'test', 'MAF', "N", 'SE', "type", "varbeta")
      df2 <- as.list(df2)
      
      
      # Convert data frames to lists
      df1_list <- list(
        snp = df1$snp,
        position = df1$position,
        chr = df1$chr,
        beta = df1$beta,
        pvalues = df1$pvalues,
        ensembl = df1$ensembl,
        Gene.name = df1$Gene.name,
        test = df1$test,
        MAF = df1$MAF,
        N = df1$N,
        type = unique(df1$type)[1]  # Ensure type is a single value
      )
      
      df2_list <- list(
        snp = df2$snp,
        position = df2$position,
        chr = df2$chr,
        beta = df2$beta,
        pvalues = df2$pvalues,
        ensembl = df2$ensembl,
        Gene.name = df2$Gene.name,
        test = df2$test,
        MAF = df2$MAF,
        N = df2$N,
        SE = df2$SE,
        type = unique(df2$type)[1],  # Ensure type is a single value
        varbeta = df2$varbeta
      )
      
      # Run coloc
      # my.res <- coloc.abf(dataset1 = df1, dataset2 = df2)
      
      my.res <- coloc.abf(dataset1 = df1_list, dataset2 = df2_list)
      
      # Get results
      results_subset <- subset(my.res$results)
      colnames(results_subset)[colnames(results_subset) == "snp"] <- "SNP"
      selected_snps <- results_subset$SNP
      selected_GWAS_eQTL <- GWAS_eQTL_trait[GWAS_eQTL_trait$SNP %in% selected_snps, ]
      
      # Merge results
      coloc_results <- merge(results_subset, selected_GWAS_eQTL, by = "SNP", all = TRUE)
      
      # Append to the trait_results data frame
      trait_results <- rbind(trait_results, coloc_results)
      
    }
      
    cat("Processed:", gene, ", ", trait_name, ", Tissue:", tissue, "\n")
    
  }
    
  # Save results for the current trait to a CSV file
  write.csv(trait_results, file = paste0('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/STARNET_results/STARNET_coloc_results_', trait_name, '.csv'), row.names = FALSE)
    
  cat("Processed trait:", trait_name, "\n")
    
  # Make sure to return the trait_results dataframe
  return(trait_results)
}

# Select phenotypes
unique_traits <- c("Body_mass_index", "Waist-hip_ratio_adj_BMI", "Waist-hip_ratio",
                   "LDL_cholesterol", "HDL_cholesterol", "Total_cholesterol", "Non-HDL_cholesterol", "Triglycerides",
                   "Fasting_insulin", "Two_hour_glucose", "Fasting_glucose", "HbA1c",
                   "Diastolic_blood_pressure", "Pulse_pressure", "Systolic_blood_pressure", 
                   "Neutrophil_count", "White_blood_cell_count", "Lymphocyte_count", "Basophil_count", "Monocyte_count", "Platelet_count", "Red_blood_cell_count", "Eosinophil_count",
                   "Type_2_diabetes", "Atrial_fibrillation", "Heart_Failure", "pad_eversmoker", "pad_nodiabetes", "pad_primary",  "NAFLD", 'Resting_heart_rate', "CAD")

# Parallel calculation
as.integer(system("nproc --all", intern = TRUE))

# Set up parallel backend with desired number of cores
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Iterate over unique trait values in parallel
# Iterate over unique trait values in parallel
trait_results_list <- foreach(trait_name = unique_traits, .combine = rbind) %dopar% {
    process_trait(trait_name)
}

# Close the parallel backend
stopCluster(cl)

# 3. Merge the coloc results
# make a subset based on HPP4 >= 0.6
# SNP.PP.H4 column
# Set your folder path
folder_path <- "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/STARNET_results/"

# Get a list of CSV files in the folder
# GTEx
csv_files <- list.files(path = folder_path, pattern = '^STARNET', full.names = TRUE)

# Initialize an empty data frame to store the merged data
merged_data <- data.frame() 
  
# Loop through the CSV files, read them, and row-bind to the merged_data
for (csv_file in csv_files) {
  data <- read.csv(csv_file, header = TRUE)
  data <- subset(data, SNP.PP.H4 >= 0.6)
  merged_data <- rbind(merged_data, data)
  
  # Print the file being processed
  cat("Processed file:", csv_file, "\n")
}  

# View the first few rows of the merged data
head(merged_data) 

write.csv(merged_data, file = "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_STARNET_eQTLs_HPP4_06_P_5e_8.csv", row.names = FALSE)

