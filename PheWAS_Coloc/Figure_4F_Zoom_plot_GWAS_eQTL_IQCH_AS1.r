library(dplyr)
library(ggplot2)
library(locuscomparer) 

# I. Prepare files

# 1. eQTL
# Select all IQCH-AS1 adipose tissue eQTLs 
GTEx_eQTLs <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/hg19_hg38_rsID_All_eQTL_GTEx_Genes_JointGenes_SNEEP_GATESGeneBodies.csv')
STARNET_eQTLs <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/All_eQTL_STARNET_JointGenes_SNEEP_GATESGeneBodies.csv')

# Make a subset by gene and tissue of interest
IQCH_AS1_GTEX_eQTLs <- subset(GTEx_eQTLs, 
                          gene_ensembl == "ENSG00000259673" & 
                          tissue %in% c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum"))

IQCH_AS1_STARNET_eQTLs <- subset(STARNET_eQTLs, 
                             gene_ensembl == "ENSG00000259673" & 
                             tissue %in% c("SF", "VAF"))

# Select the columns with rsID and p-value for Zoom plot
IQCH_AS1_GTEX_eQTLs_ZP <- IQCH_AS1_GTEX_eQTLs %>%
  select(rsid = rsID, pval = pval_nominal, position = position_hg19) %>%
  distinct()

IQCH_AS1_STARNET_eQTLs_ZP <- IQCH_AS1_STARNET_eQTLs %>%
  select(rsid = SNP, pval = p.value, position = pos) %>%
  distinct()

IQCH_AS1_all_adipose_eQTL_ZP <- rbind(IQCH_AS1_GTEX_eQTLs_ZP, IQCH_AS1_STARNET_eQTLs_ZP)

# Make a subset for IQCH-AS1
# Position: hg19 chr15:67,807,260-67,813,400 
IQCH_AS1_all_adipose_eQTL_ZP <- subset(IQCH_AS1_all_adipose_eQTL_ZP, position >= 67500000 & position <= 68200000 & chr == 15)

IQCH_AS1_all_adipose_eQTL_ZP$rsid <- as.character(IQCH_AS1_all_adipose_eQTL_ZP)
IQCH_AS1_all_adipose_eQTL_ZP$position <- NULL

write.table(IQCH_AS1_all_adipose_eQTL_ZP, file = "Zoom_Plot_IQCH_AS1/Adipose_eQTLs_IQCH_AS1_locuscomparer.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 2. GWAS
# Select SNPs within the IQCH-AS1 region from the WHR GWAS 
whr_trait <- read.table("bmi_whr/giant-ukbb.meta-analysis.combined.23May2018__Waist-hip_ratio.txt.gz", header = TRUE, stringsAsFactors = FALSE)

# Make a subset for IQCH-AS1
# Position: hg19 chr15:67,807,260-67,813,400 
IQCH_AS1_WHR_ZP <- subset(whr_trait, position >= 67500000 & position <= 68200000 & chr == 15)

IQCH_AS1_WHR_ZP <- whr_trait %>%
  select(rsid = SNP, pval = P) %>%
  distinct()

# Remove parts after the first column
IQCH_AS1_WHR_ZP$rsid <- sub(":.*", "", IQCH_AS1_WHR_ZP$rsid)
IQCH_AS1_WHR_ZP$rsid <- as.character(IQCH_AS1_WHR_ZP$rsid)

write.table(IQCH_AS1_WHR_ZP, file = "Zoom_Plot_IQCH_AS1/WHR_IQCH_AS1_locuscomparer.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# II. Locus comparer zoom plots
# Lead WHR SNP
gwas_fn = 'Zoom_Plot_IQCH_AS1/WHR_IQCH_AS1_locuscomparer.tsv'
eqtl_fn = 'Zoom_Plot_IQCH_AS1/Adipose_eQTLs_IQCH_AS1_locuscomparer.tsv'

LZ_plot <- locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'WHR GWAS', title2 = 'Adipose eQTL', snp = 'rs3784699')

ggsave(filename = 'Plots/Figure_4/Figure_4F/WHR_eQTL_locuscomparer_plot_rs3784699.png', plot = LZ_plot, width = 8, height = 3.8, dpi = 800)

pdf_file <- 'Plots/Figure_4/Figure_4F/WHR_eQTL_locuscomparer_plot_rs3784699.pdf'
ggsave(filename = pdf_file, plot = LZ_plot, width = 8, height = 3.8, dpi = 800)

LZ_plot

