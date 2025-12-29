library(dplyr)
library(ggplot2)
library(locuscomparer) 

# I. Prepare files

# 1. eQTL
# Select VAF eQTLs from STARNET
STARNET_eQTLs <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/All_eQTL_STARNET_JointGenes_SNEEP_GATESGeneBodies.csv')

# Make a subset by gene and tissue of interest
IQCH_STARNET_eQTLs <- subset(STARNET_VAF_eQTLs, gene_ensembl == "ENSG00000259673")
IQCH_STARNET_eQTLs <- subset(IQCH_STARNET_eQTLs, tissue == "VAF")

head(IQCH_STARNET_eQTLs)

# Select the column with rsID and p-value
IQCH_AS1_VAF_eQTLs_ZP_STARNET <- IQCH_STARNET_eQTLs[, c("SNP", "p.value")]
colnames(IQCH_AS1_VAF_eQTLs_ZP_STARNET) <- c("rsid", "pval")
IQCH_AS1_VAF_eQTLs_ZP_STARNET <- IQCH_AS1_VAF_eQTLs_ZP_STARNET %>% distinct()
IQCH_AS1_VAF_eQTLs_ZP_STARNET <- subset(IQCH_AS1_VAF_eQTLs_ZP_STARNET, rsid != 'MERGED_DEL_2_80211')

IQCH_AS1_VAF_eQTLs_ZP_STARNET$rsid <- as.character(IQCH_AS1_VAF_eQTLs_ZP_STARNET$rsid)
head(IQCH_AS1_VAF_eQTLs_ZP_STARNET)

write.table(IQCH_AS1_VAF_eQTLs_ZP_STARNET, file = "Zoom_Plot_IQCH_AS1/VAF_STARNET_eQTLs_IQCH_AS1_locuscomparer.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# 2. GWAS
# Select SNPs within the IQCH-AS1 region from the BMI GWAS 
bmi_trait <- read.table("bmi_whr/giant-ukbb.meta-analysis.combined.23May2018__Body_mass_index.txt.gz", header = TRUE, stringsAsFactors = FALSE)

# Make a subset for IQCH-AS1
# Position: hg19 chr15:67,807,260-67,813,400 
IQCH_AS1_BMI_ZP <- subset(bmi_trait, position >= 67500000 & position <= 68200000 & chr == 15)

IQCH_AS1_BMI_ZP <- bmi_trait %>%
  select(rsid = SNP, pval = P) %>%
  distinct()

# Remove parts after the first column
IQCH_AS1_BMI_ZP$rsid <- sub(":.*", "", IQCH_AS1_BMI_ZP$rsid)
IQCH_AS1_BMI_ZP$rsid <- as.character(IQCH_AS1_BMI_ZP$rsid)

write.table(IQCH_AS1_BMI_ZP, file = "Zoom_Plot_IQCH_AS1/BMI_IQCH_AS1_locuscomparer.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# II. Locus comparer zoom plots
# rs4261477 - lead BMI SNP
gwas_fn = 'Zoom_Plot_IQCH_AS1/BMI_IQCH_AS1_locuscomparer.tsv'
eqtl_fn = 'Zoom_Plot_IQCH_AS1/VAF_STARNET_eQTLs_IQCH_AS1_locuscomparer.tsv'

LZ_plot <-locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'BMI GWAS', title2 = 'VAF eQTL', snp = 'rs4261477')

ggsave(filename = 'Plots/Figure_4/Figure_4F/BMI_VAF_eQTL_locuscompare_plot_rs4261477.png', plot = LZ_plot, width = 7, height = 4, dpi = 800)

pdf_file <- 'Plots/Figure_4/Figure_4F/BMI_VAF_eQTL_locuscompare_plot_rs4261477.pdf'
ggsave(filename = pdf_file, plot = LZ_plot, width = 7, height = 4, dpi = 800)

LZ_plot
