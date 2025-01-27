library(tidyverse)
library(data.table)
library(stringr)
library(rtracklayer)

# Read eQTL processed files
# STARNET 
STARNET_eQTL <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/All_eQTL_STARNET_JointGenes_SNEEP_GATESGeneBodies.csv')
STARNET_eQTL$SNP <- NULL
names(STARNET_eQTL)[names(STARNET_eQTL) == "snpid"] <- "SNP"
STARNET_eQTL$chr_bp_position_hg19 <- paste(STARNET_eQTL$chr, STARNET_eQTL$pos, sep = ":")
STARNET_eQTL$markername <- NULL

# GTEx
GTEx_eQTL <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/hg19_hg38_rsID_All_eQTL_GTEx_Genes_JointGenes_SNEEP_GATESGeneBodies.csv')

# Create a list of data frames
data_frames <- list(GTEx_eQTL, STARNET_eQTL)

output_path <- 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/'

# I.
# GWAS data extraction

# Heart failure
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/GWAS_Catalog/HF/GCST90162626_buildGRCh37.tsv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

trait_name <- 'Heart_Failure'
rs_id_colname <- 'variant_id'

# Iterate through the list
for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")

  # Set eQTL_file to the current data frame
  eQTL_file <- df
  
  # Find shared SNPs
  shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
  
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
  
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "\n")
  
}

# Lipids
# List of trait names and corresponding file names
trait_file_pairs <- list(
    c("HDL_cholesterol", '/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_lipids/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_HDL_INV_ALL_with_N_1__HDL_cholesterol.gz'),
  c("LDL_cholesterol", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_lipids/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_LDL_INV_ALL_with_N_1__LDL_cholesterol.gz"),
  c("Triglycerides", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_lipids/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_logTG_INV_ALL_with_N_1__Triglycerides.gz"),
  c("Non-HDL_cholesterol", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_lipids/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_nonHDL_INV_ALL_with_N_1__Non-HDL_cholesterol.gz"),
  c("Total_cholesterol", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_lipids/with_BF_meta-analysis_AFR_EAS_EUR_HIS_SAS_TC_INV_ALL_with_N_1__Total_cholesterol.gz")
)

rs_id_colname <- 'rsID'

# Iterate through trait names and file names
for (trait_file_pair in trait_file_pairs) {
  trait_name <- trait_file_pair[1]
  gwas_file <- trait_file_pair[2]
  
  # Read the GWAS data for the current trait
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  
  for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
  
    # Set eQTL_file to the current data frame
    eQTL_file <- df
  
    # Find shared SNPs
    shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "\n")
  }
  
}

# Glucose
# List of trait names and corresponding file names
trait_file_pairs <- list(
  c("Two_hour_glucose", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/MAGIC1000G_2hGlu_EUR.tsv.gz"),
  c("Fasting_glucose", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/MAGIC1000G_FG_EUR.tsv.gz"),
  c("Fasting_insulin", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/MAGIC1000G_FI_EUR.tsv.gz"),
  c("HbA1c", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/MAGIC1000G_HbA1c_EUR.tsv.gz")
 )

rs_id_colname <- 'variant'

# Iterate through trait names and file names
for (trait_file_pair in trait_file_pairs) {
    trait_name <- trait_file_pair[1]
    gwas_file <- trait_file_pair[2]
    
    # Read the GWAS data for the current trait
    gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
    
    for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
  
    # Set eQTL_file to the current data frame
    eQTL_file <- df
    
    # Find shared SNPs
    shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "\n")
  }
  
}

# Anthropometric traits
# List of trait names and corresponding file names
trait_file_pairs <- list(
  c("Waist-hip_ratio_adj_BMI", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/giant-ukbb.meta-analysis.combined.23May2018__Waist-hip_ratio_adj_BMI.txt.gz"),
  
    c("Body_mass_index", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/bmi_whr/giant-ukbb.meta-analysis.combined.23May2018__Body_mass_index.txt.gz"),
  
    c("Waist-hip_ratio", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/bmi_whr/giant-ukbb.meta-analysis.combined.23May2018__Waist-hip_ratio.txt.gz")
)

rs_id_colname <- 'RSID'

# Iterate through trait names and file names
for (trait_file_pair in trait_file_pairs) {
  trait_name <- trait_file_pair[1]
  gwas_file <- trait_file_pair[2]
  
  # Read the GWAS data for the current trait
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  gwas$RSID <- unlist(lapply(strsplit(gwas$SNP, ':'), FUN=function(x){x[1]}))
  
  for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
  
    # Set eQTL_file to the current data frame
    eQTL_file <- df
  
    # Find shared SNPs
    shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "\n")
  }
  
}

# T2 Diabetes
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/diabetes/DIAMANTE-TA_T2D_1.2million.sumstat__Type_2_diabetes.txt.gz", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

#head(gwas)
trait_name <- 'Type_2_diabetes'
rs_id_colname <- 'rsID'

  for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")

  # Set eQTL_file to the current data frame
  eQTL_file <- df
  
  # Find shared SNPs
  shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
  
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
  
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "\n")
  
}

# Nonalcoholic fatty liver disease (NAFLD)
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/GCST90091033_buildGRCh37__Non-alcoholic_fatty_liver_disease.tsv.gz", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

# head(gwas)
trait_name <- 'NAFLD'
rs_id_colname <- 'variant_id'

  for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")

  # Set eQTL_file to the current data frame
  eQTL_file <- df
  
  # Find shared SNPs
  shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
  
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
  
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "\n")
  
}

# Atrial fibrillation
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/AF_HRC_GWAS_ALLv11_29892015__0.6million__Atrial_fibrillation.txt", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

trait_name <- 'Atrial_fibrillation'
rs_id_colname <- 'MarkerName'

  for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")

  # Set eQTL_file to the current data frame
  eQTL_file <- df
  
  # Find shared SNPs
  shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
  
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
  
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "\n")
  
}

# Peripheral artery disease
# Get a list of files in the specified directory
file_dir <- "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/vanzuydametal.2020_PAD/"
file_paths <- list.files(file_dir, pattern = "\\.txt\\.gz", full.names = TRUE)

rs_id_colname <- 'snp'

# Iterate through the list of file paths
for (gwas_file in file_paths) {
  # Extract trait_name from the filename (remove directory path and extension)
  trait_name <- sub(".txt.gz$", "", basename(gwas_file))
  
  # Read the GWAS data for the current trait
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  
    for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
  
    # Set eQTL_file to the current data frame
    eQTL_file <- df  
    
    # Find shared SNPs
    shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "\n")
  }
}

# Heart rate
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc__Resting_heart_rate.gz", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
gwas$sample_size <- NA

trait_name <- 'Resting_heart_rate'
rs_id_colname <- 'SNP_GWAS'

for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
    
  # Set eQTL_file to the current data frame
  eQTL_file <- df    

  # Find shared SNPs
  shared_snps <- intersect(gwas[[rs_id_colname]], eQTL_file$SNP)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas[[rs_id_colname]] %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file[eQTL_file$SNP %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = rs_id_colname, by.y = "SNP")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
  
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
  
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "\n")
  
} 

# Blood cells 
# List of trait names and corresponding file names
trait_file_pairs <- list(  
  c("Basophil_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_BAS_Trans_GWAMA.out__Basophil_count.gz"),
  c("Eosinophil_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_EOS_Trans_GWAMA.out__Eosinophil_count.gz"),
  c("Lymphocyte_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_LYM_Trans_GWAMA.out__Lymphocyte_count.gz"),
  c("Monocyte_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_MON_Trans_GWAMA.out__Monocyte_count.gz"),
  c("Neutrophil_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_NEU_Trans_GWAMA.out__Neutrophil_count.gz"),
  c("Platelet_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_PLT_Trans_GWAMA.out__Platelet_count.gz"),
  c("Red_blood_cell_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_RBC_Trans_GWAMA.out__Red_blood_cell_count.gz"),
  c("White_blood_cell_count", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/blood_immune_cells/BCX2_WBC_Trans_GWAMA.out__White_blood_cell_count.gz")
)


# Iterate through trait names and file names
for (trait_file_pair in trait_file_pairs) {
  trait_name <- trait_file_pair[1]
  gwas_file <- trait_file_pair[2]
  
  # Read the GWAS data for the current trait
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  gwas$chr_bp_position <- sub("_.*", "", gwas$rs_number)

  # Iterate through the list of data frames
  for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
    
    # Set eQTL_file to the current data frame
    eQTL_file_hg19 <- df    
    
    # Find shared SNPs
    shared_snps <- intersect(gwas$chr_bp_position, eQTL_file_hg19$chr_bp_position_hg19)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas$chr_bp_position %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file_hg19[eQTL_file_hg19$chr_bp_position_hg19 %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = "chr_bp_position", by.y = "chr_bp_position_hg19")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "for trait:", trait_name, "\n")
  }
}

# Blood pressure
# List of trait names and corresponding file names
trait_file_pairs <- list(
  c("Diastolic_blood_pressure", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/Evangelou_30224653_1million__Diastolic_blood_pressure.txt.gz"),
  c("Pulse_pressure", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/Evangelou_30224653_1million__Pulse_pressure.txt.gz"),
  c("Systolic_blood_pressure", "/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/Evangelou_30224653_1million__Systolic_blood_pressure.txt.gz")
)

# Iterate through trait names and file names
for (trait_file_pair in trait_file_pairs) {
  trait_name <- trait_file_pair[1]
  gwas_file <- trait_file_pair[2]
  
  # Read the GWAS data for the current trait
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
  gwas$chr_bp_position <- sub("^(.*?:[^:]+):.*", "\\1", gwas$MarkerName)
  names(gwas)[names(gwas) == "P"] <- "p_value_GWAS"
  names(gwas)[names(gwas) == "Effect"] <- "beta_gwas"
  names(gwas)[names(gwas) == "N_effective"] <- "sample_size"
  
  # Iterate through the list of data frames
  for (df in data_frames) {
    # Define the name based on the current data frame
    df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
    
    # Set eQTL_file to the current data frame
    eQTL_file_hg19 <- df      
    # Find shared SNPs
    shared_snps <- intersect(gwas$chr_bp_position, eQTL_file_hg19$chr_bp_position_hg19)
    
    # Extract shared SNPs
    gwas_subset <- gwas[gwas$chr_bp_position %in% shared_snps, ]
    eQTL_file_subset <- eQTL_file_hg19[eQTL_file_hg19$chr_bp_position_hg19 %in% shared_snps, ]
    
    # Merge subsets and add Trait column
    merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = "chr_bp_position", by.y = "chr_bp_position_hg19")
    merged_subsets$Trait <- trait_name
    
    # Define the output file name
    output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
    # Save the merged data to a CSV file
    write.csv(merged_subsets, output_file)
    
    # Print a message to indicate completion
    cat("Processed data frame:", df_name, "for trait:", trait_name, "\n")
  }
}

# Heart rate
# Read GWAS sum stat
gwas <- read.table("/storage2/Public_data/PheWAS/39_CAD_relevant_summary_statisitcis/original/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc__Resting_heart_rate.gz", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

gwas$chr_bp_position <- sub("_.*", "", gwas$SNP)
trait_name <- 'Resting_heart_rate_bp'

# Iterate through the list of data frames
for (df in data_frames) {
  # Define the name based on the current data frame
  df_name <- ifelse(identical(df, GTEx_eQTL), "GTEx", "STARNET")
    
  # Set eQTL_file to the current data frame
  eQTL_file_hg19 <- df    
    
  # Find shared SNPs
  shared_snps <- intersect(gwas$chr_bp_position, eQTL_file_hg19$chr_bp_position_hg19)
  
  # Extract shared SNPs
  gwas_subset <- gwas[gwas$chr_bp_position %in% shared_snps, ]
  eQTL_file_subset <- eQTL_file_hg19[eQTL_file_hg19$chr_bp_position_hg19 %in% shared_snps, ]
  
  # Merge subsets and add Trait column
  merged_subsets <- merge(gwas_subset, eQTL_file_subset, by.x = "chr_bp_position", by.y = "chr_bp_position_hg19")
  merged_subsets$Trait <- trait_name
  
  # Define the output file name
  output_file <- paste0(output_path, df_name, '_', trait_name, '.csv')
    
  # Save the merged data to a CSV file
  write.csv(merged_subsets, output_file)
    
  # Print a message to indicate completion
  cat("Processed data frame:", df_name, "for trait:", trait_name, "\n")
}


# II
# Get uni columns and merge all the files together
file_path <- 'results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/'

# Atrial_fibrillation
# GTEx then STARNET
trait_name <- 'Atrial_fibrillation'

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_name), paste0("GTEx_", trait_name))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "MarkerName"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "P.value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "Effect"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "chr.x"] <- "chr"
    
    file_by_trait$sample_size <- 588190
    file_by_trait$freq <- NA
    names(file_by_trait)[names(file_by_trait) == "StdErr"] <- "SE"


    # Make a subset of unified columns for GTEx
   file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
      
    } else {
      # Rename columns for STARNET
      names(file_by_trait)[names(file_by_trait) == "MarkerName"] <- "SNP"
      names(file_by_trait)[names(file_by_trait) == "P.value"] <- "p_value_GWAS"
      names(file_by_trait)[names(file_by_trait) == "Effect"] <- "beta_GWAS"
      names(file_by_trait)[names(file_by_trait) == "chr.x"] <- "chr"
      names(file_by_trait)[names(file_by_trait) == "pos.x"] <- "position_hg19"
      names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
      names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
      names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
      names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"

      file_by_trait$sample_size <- 588190
      file_by_trait$freq <- NA
      names(file_by_trait)[names(file_by_trait) == "StdErr"] <- "SE"
      
      file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}


# Anthropometrics traits
trait_names <- c(
  "Waist-hip_ratio_adj_BMI", "Body_mass_index", "Waist-hip_ratio")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
  # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "RSID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "P"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "BETA"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "N"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "Tested_Allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "Other_Allele"] <- "Allele2"
    
    names(file_by_trait)[names(file_by_trait) == "CHR"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "Freq_Tested_Allele"] <- "freq"   
    
    # Make a subset of unified columns for GTEx
   file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET
    names(file_by_trait)[names(file_by_trait) == "RSID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "P"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "BETA"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "N"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "Tested_Allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "Other_Allele"] <- "Allele2"
    
    names(file_by_trait)[names(file_by_trait) == "CHR"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "Freq_Tested_Allele"] <- "freq"
    
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL" 
      
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Blood
trait_names <- c("Basophil_count", "Eosinophil_count", "Lymphocyte_count", "Monocyte_count", "Neutrophil_count", "Platelet_count", "Red_blood_cell_count", "White_blood_cell_count")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "RSID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "n_samples"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "reference_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    names(file_by_trait)[names(file_by_trait) == "chr_bp_position"] <- "chr_bp_position_hg19"
 
    names(file_by_trait)[names(file_by_trait) == "se"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "chr_name_hg19"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "eaf"] <- "freq"

    # Make a subset of unified columns for GTEx
    file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET    
    names(file_by_trait)[names(file_by_trait) == "RSID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p.value.x"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta.x"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "n_samples"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "reference_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    names(file_by_trait)[names(file_by_trait) == "chr_bp_position"] <- "chr_bp_position_hg19"
 
    names(file_by_trait)[names(file_by_trait) == "se"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "chr_name_hg19"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "eaf"] <- "freq"

    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta.y"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value.y"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"   
    
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}


# Peripheral artery disease
trait_names <- c('pad_diabetes', "pad_neversmoker", "pad_eversmoker", "pad_nodiabetes", "pad_primary")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "snp"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p"] <- "p_value_GWAS"
    
    file_by_trait <- file_by_trait %>% mutate(beta_GWAS = log(or, base = exp(1)))
    
    file_by_trait$sample_size <- file_by_trait$ncases + file_by_trait$ncontrols
    names(file_by_trait)[names(file_by_trait) == "or_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "nonor_allele"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "or_se"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "eaf"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "chr.x"] <- "chr"

    # Make a subset of unified columns for GTEx
    file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET    
    names(file_by_trait)[names(file_by_trait) == "snp"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p"] <- "p_value_GWAS"
    
    file_by_trait <- file_by_trait %>% mutate(beta_GWAS = log(or, base = exp(1)))
    
    file_by_trait$sample_size <- file_by_trait$ncases + file_by_trait$ncontrols
    names(file_by_trait)[names(file_by_trait) == "or_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "nonor_allele"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "or_se"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "eaf"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "chr.x"] <- "chr"
    
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"     
    
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Lipids
trait_names <- c("HDL_cholesterol", "LDL_cholesterol", "Triglycerides", "Non-HDL_cholesterol", "Total_cholesterol")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "pvalue"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "or"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "N"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "ALT"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "REF"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "METAL_StdErr"] <- "SE"
    
    names(file_by_trait)[names(file_by_trait) == "CHROM"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "POOLED_ALT_AF"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "METAL_Effect"] <- "beta_GWAS"
    
    # Make a subset of unified columns for GTEx
    file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET  
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "pvalue"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "or"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "N"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "ALT"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "REF"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "METAL_StdErr"] <- "SE"
    
    names(file_by_trait)[names(file_by_trait) == "CHROM"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "POOLED_ALT_AF"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "METAL_Effect"] <- "beta_GWAS"
    
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"     
    
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Glucose
trait_names <- c("Two_hour_glucose", "Fasting_glucose", "Fasting_insulin", "HbA1c")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "variant"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"
    
    # Make a subset of unified columns for GTEx
    file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET  
    names(file_by_trait)[names(file_by_trait) == "variant"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta.x"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"

    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"
    
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta.y"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"  
    
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Pressure
trait_names <- c("Diastolic_blood_pressure", "Pulse_pressure", "Systolic_blood_pressure")

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_names), paste0("GTEx_", trait_names))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "beta_gwas"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "StdErr"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "Freq1"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "chr_name_hg19"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "chr_bp_position"] <- "chr_bp_position_hg19"
    
    # Make a subset of unified columns for GTEx
    file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
   
    } else {
    # Rename columns for STARNET   
    names(file_by_trait)[names(file_by_trait) == "beta_gwas"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "StdErr"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "Freq1"] <- "freq"
    names(file_by_trait)[names(file_by_trait) == "chr_name_hg19"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "chr_bp_position"] <- "chr_bp_position_hg19"
    
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL"     
    
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# T2D
trait_name <- 'Type_2_diabetes'
# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_name), paste0("GTEx_", trait_name))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_p.value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_beta"] <- "beta_GWAS"
    
    names(file_by_trait)[names(file_by_trait) == "Effective_sample_size"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    
    names(file_by_trait)[names(file_by_trait) == "chromosome.b37."] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_SE"] <- "SE"
    
    file_by_trait$freq <- NA
    
    # Make a subset of unified columns for GTEx
  file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
      
    } else {
    # Rename columns for STARNET
    names(file_by_trait)[names(file_by_trait) == "rsID"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_p.value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_beta"] <- "beta_GWAS"
    
    names(file_by_trait)[names(file_by_trait) == "Effective_sample_size"] <- "sample_size"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    
    names(file_by_trait)[names(file_by_trait) == "chromosome.b37."] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "Fixed.effects_SE"] <- "SE"  
      
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL" 
    
    file_by_trait$freq <- NA
      
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Heart_Failure
trait_name <- 'Heart_Failure'
# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_name), paste0("GTEx_", trait_name))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    names(file_by_trait)[names(file_by_trait) == "variant_id"] <- "SNP"
    
    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"
    
    names(file_by_trait)[names(file_by_trait) == "variant_id.y"] <- "variant_id"
    
    file_by_trait$sample_size <- 1665481

    # Make a subset of unified columns for GTEx
   file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
      
    } else {
    # Rename columns for STARNET
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta.x"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    names(file_by_trait)[names(file_by_trait) == "variant_id"] <- "SNP"
    
    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"
    
    file_by_trait$sample_size <- 1665481
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta.y"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL" 
      
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Nonalcoholic fatty liver disease (NAFLD)
trait_name <- 'NAFLD'
# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_name), paste0("GTEx_", trait_name))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "variant_id"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    file_by_trait$sample_size <- 778614
    
    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"
    
    names(file_by_trait)[names(file_by_trait) == "variant_id.y"] <- "variant_id"

    # Make a subset of unified columns for GTEx
   file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
      
    } else {
    # Rename columns for STARNET
     names(file_by_trait)[names(file_by_trait) == "variant_id"] <- "SNP"
    names(file_by_trait)[names(file_by_trait) == "p_value"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "beta.x"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "effect_allele"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "other_allele"] <- "Allele2"
    file_by_trait$sample_size <- 778614
    
    names(file_by_trait)[names(file_by_trait) == "chromosome"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "standard_error"] <- "SE"
    names(file_by_trait)[names(file_by_trait) == "effect_allele_frequency"] <- "freq"     
      
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta.y"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL" 
      
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}

# Resting heart rate
trait_name <- 'Resting_heart_rate'

# Add prefixes 'GTEx_' and 'STARNET_'
prefixed_trait_names <- c(paste0("STARNET_", trait_name), paste0("GTEx_", trait_name))

# Iterate through the list of file name prefixes
for (prefix in prefixed_trait_names) {
  # Create file name and full file path
  full_file_path <- file.path(file_path, paste0(prefix, '.csv'))

  # Read the CSV file
  file_by_trait <- read.csv(full_file_path)

  if (startsWith(prefix, "GTEx_")) {
    # Rename columns for GTEx
    names(file_by_trait)[names(file_by_trait) == "P"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "BETA"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "A1"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "A0"] <- "Allele2"    
    names(file_by_trait)[names(file_by_trait) == "CHR"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "MAF"] <- "freq"
  
    file_by_trait$sample_size <- 458969

    # Make a subset of unified columns for GTEx
   file_by_trait <- file_by_trait[, c("gene_ensembl", "SNP", "gene_id", "chr", "chr_bp_position_hg19", 'chr_pos_hg38', "slope", "slope_se", 'ma_samples', 'ma_count', 'pval_nominal', "pval_beta", "pval_nominal_threshold", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2", 'freq',"SE", 'sample_size', "tissue","Trait", "test", "Gene.name", "Biotype")]
      
    } else {
    # Rename columns for STARNET
    names(file_by_trait)[names(file_by_trait) == "P"] <- "p_value_GWAS"
    names(file_by_trait)[names(file_by_trait) == "BETA"] <- "beta_GWAS"
    names(file_by_trait)[names(file_by_trait) == "A1"] <- "Allele1"
    names(file_by_trait)[names(file_by_trait) == "A0"] <- "Allele2"   
    names(file_by_trait)[names(file_by_trait) == "CHR"] <- "chr"
    names(file_by_trait)[names(file_by_trait) == "MAF"] <- "freq"  
    file_by_trait$sample_size <- 458969
      
    names(file_by_trait)[names(file_by_trait) == "pos"] <- "position_hg19"
    names(file_by_trait)[names(file_by_trait) == "beta"] <- "beta_eQTL"
    names(file_by_trait)[names(file_by_trait) == "p.value"] <- "p_value_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_fdr"] <- "padj_fdr_eQTL"
    names(file_by_trait)[names(file_by_trait) == "padj_bonferroni"] <- "padj_bonferroni_eQTL" 
      
    file_by_trait <- file_by_trait[, c('gene_ensembl', "SNP", "chr", "position_hg19", "chr_bp_position_hg19", "beta_eQTL", "t.stat", "p_value_eQTL", "padj_fdr_eQTL", "padj_bonferroni_eQTL", "p_value_GWAS", "beta_GWAS", "Allele1", "Allele2",  'freq', "SE", 'sample_size', "tissue", "Trait", "test", "Gene.name", "Biotype")]
    }     


  # Define the output path and write the CSV
  output_path <- file.path(file_path, 'Subsetted_df_uni_columns', paste0(prefix, '_uni_columns.csv'))
  write.csv(file_by_trait, file = output_path, row.names = FALSE)
  
  cat("Processed file:", prefix, "\n")
}


# III 
# Merge all GTEx and all STARNET subsetted files with uni colnames for heatmaps

# 1
# For STARNET
# 5e-5
# Set your folder path
folder_path <- "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Subsetted_df_uni_columns"

# Get a list of CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = '^STARNET', full.names = TRUE)

# Initialize an empty data frame to store the merged data
merged_data <- data.frame()

# Loop through the CSV files, read them, and row-bind to the merged_data
for (csv_file in csv_files) {
  data <- read.csv(csv_file, header = TRUE)
  data <- subset(data, p_value_GWAS < 5e-5)
  merged_data <- rbind(merged_data, data)
  
  # Print the file being processed
  cat("Processed file:", csv_file, "\n")
}

write.csv(merged_data, file = "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Final_merged_tables/PheWAS_results_STARNET_eQTLs_FINAL_Genes_Joint_GeneBase_5e_5.csv", row.names = FALSE)


# 2
# For GTEx
# 5e-5
# Set your folder path
folder_path <- "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Subsetted_df_uni_columns"

# Get a list of CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = '^GTEx', full.names = TRUE)

# Initialize an empty data frame to store the merged data
merged_data <- data.frame() 
 
# Loop through the CSV files, read them, and row-bind to the merged_data
for (csv_file in csv_files) {
  data <- read.csv(csv_file, header = TRUE)
  data <- subset(data, p_value_GWAS < 5e-5)
  merged_data <- rbind(merged_data, data)
  
  # Print the file being processed
  cat("Processed file:", csv_file, "\n")
}  

write.csv(merged_data, file = "results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Final_merged_tables/PheWAS_results_GTEx_eQTLs_FINAL_Genes_Joint_GeneBase_5e_5.csv", row.names = FALSE)

