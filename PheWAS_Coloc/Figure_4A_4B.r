library(dplyr)  
library(readxl)
library(ggplot2)
library(plyr)
library(Cairo)

coloc_output <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_STARNET_merged_eQTLs_HPP4_06_P_5e_8_with_eQTL_alleles_Renamed_pheno_renamed_tissue.csv')

# Perform grouping by traits
coloc_output$Trait_merged_groups <- coloc_output$Trait

coloc_output$Trait_merged_groups <- revalue(coloc_output$Trait_merged_groups,
                                   c("Atrial_fibrillation" = "Atrial fibrillation",
                                     "Heart_Failure" = "Heart failure",
                                     'Resting_heart_rate' = 'Resting heart rate',
                                     "Type_2_diabetes" = "Glycemic traits",
                                     "NAFLD" = "Nonalcoholic fatty liver disease",
                                     "HDL_cholesterol" = "Lipids",
                                     "LDL_cholesterol" = "Lipids",
                                     "Triglycerides" = "Lipids",
                                     "Total_cholesterol" = "Lipids",
                                     "Non-HDL_cholesterol" = "Lipids",
                                     "Fasting_glucose" = "Glycemic traits",
                                     "Fasting_insulin" = "Glycemic traits",
                                     "HbA1c" = "Glycemic traits",
                                     "Diastolic_blood_pressure" = "Blood pressure",
                                     "Systolic_blood_pressure" = "Blood pressure",
                                     "Pulse_pressure" = "Blood pressure",
                                     "Two_hour_glucose" = "Glycemic traits",
                                     "Basophil_count" = "Inflammatory biomarkers",
                                     "Eosinophil_count" = "Inflammatory biomarkers",
                                     "Lymphocyte_count" = "Inflammatory biomarkers",
                                     "Monocyte_count" = "Inflammatory biomarkers",
                                     "Neutrophil_count" = "Inflammatory biomarkers",
                                     "Red_blood_cell_count" = "Inflammatory biomarkers",
                                     "Platelet_count" = "Coagulation/Thrombosis",
                                     "White_blood_cell_count" = "Inflammatory biomarkers",
                                     "Body_mass_index" = "Adiposity",
                                     "Waist-hip_ratio" = "Adiposity",
                                     "Waist-hip_ratio_adj_BMI" = "Adiposity",
                                     "pad_primary" = "Peripheral arterial disease",
                                     "pad_diabetes" = "Peripheral arterial disease",
                                     "pad_nodiabetes" = "Peripheral arterial disease",
                                     "pad_eversmoker" = "Peripheral arterial disease",
                                     "pad_neversmoker" = "Peripheral arterial disease"
                                     ))

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
coloc_output <- coloc_output %>%
  left_join(gene_table %>% select(Gene.name, Conserved.in.mouse), by = "Gene.name")

# 'known.CAD.loci.gene' column
coloc_output <- coloc_output %>%
  left_join(gene_table %>% select(Gene.name, known_CAD_loci_gene_new_list), by = "Gene.name")

head(coloc_output)

# I.
# Figure 4A
# All CAD candidate genes

# 4A - 1
# Gene freq by phenotype 
plot_data <- coloc_output
plot_data <- plot_data %>%
  dplyr::group_by(Trait_merged_groups, Biotype_binary) %>%
  dplyr::summarise(unique_genes = dplyr::n_distinct(Gene.name)) %>%
  dplyr::ungroup()

# Calculate the total number of unique genes per trait and sort by it
plot_data <- plot_data %>%
  group_by(Trait_merged_groups) %>%
  mutate(total_unique_genes = sum(unique_genes)) %>%
  ungroup() %>%
  arrange(desc(total_unique_genes))

trait_order <- plot_data %>%
  dplyr::group_by(Trait_merged_groups) %>%
  dplyr::summarise(total_unique_genes = sum(unique_genes)) %>%
  dplyr::arrange(desc(total_unique_genes)) %>%
  dplyr::pull(Trait_merged_groups)

# Ensure that the Trait_merged_groups levels are ordered from high to low based on the total number of unique genes
plot_data$Trait_merged_groups <- factor(plot_data$Trait_merged_groups, levels = trait_order)

# Create a barplot
plot_data$Biotype_binary <- factor(plot_data$Biotype_binary, levels = c("protein-coding", "non-coding"))

stacked_bar_plot <- ggplot(plot_data, aes(x = Trait_merged_groups, y = unique_genes, fill = Biotype_binary)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("protein-coding" = "#FFC20A", "non-coding" = "#0C7BDC")) +
  labs(title = '', x = "Phenotype", y = "Number of genes", fill = "Gene biotype") +
  scale_y_continuous(breaks = seq(0, 1000, by = 200), limits = c(0, 1010)) +  # Set y-axis limits and labels
  theme_classic() +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, family = "Calibri", color = "black"),
    axis.text.y = element_text(size = 8, family = "Calibri", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 9, family = "Calibri", face = "bold", color = "black"),  # Enlarge plot title
    axis.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge axis labels
    legend.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge legend title
    legend.text = element_text(size = 8, family = "Calibri", color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(size = 0.3)) +
  
    theme(legend.key.height = unit(0.2, "mm"), legend.key.width = unit(2, "mm"))  # Adjust size of legend keys

stacked_bar_plot

ggsave(stacked_bar_plot, filename = "Figure_4/Figure_4A/Fig_4A_Barplot_pheno_associat_All_gene_count_by_pheno.pdf", width = 3.05, height = 2.88, device = cairo_pdf, dpi = 800)

# 4A - 2
# Gene count by tissue
plot_data_tissue <- coloc_output

# Group by tissue and count
plot_data_tissue <- plot_data_tissue %>%
  dplyr::group_by(tissue_labels, Biotype_binary) %>%
  dplyr::summarise(unique_genes = dplyr::n_distinct(Gene.name)) %>%
  dplyr::ungroup()

# Calculate the total number of unique genes per tissue and sort by it
plot_data_tissue <- plot_data_tissue %>%
  group_by(tissue_labels) %>%
  mutate(total_unique_genes = sum(unique_genes)) %>%
  ungroup() %>%
  arrange(desc(total_unique_genes))

tissue_order <- plot_data_tissue %>%
  dplyr::group_by(tissue_labels) %>%
  dplyr::summarise(total_unique_genes = sum(unique_genes)) %>%
  dplyr::arrange(desc(total_unique_genes)) %>%
  dplyr::pull(tissue_labels)

# Ensure that the tissue_labels levels are ordered from high to low based on the total number of unique genes
plot_data_tissue$tissue_labels <- factor(plot_data_tissue$tissue_labels, levels = tissue_order)

# Create a barplot
# Convert Biotype_binary to a factor with the desired order
plot_data_tissue$Biotype_binary <- factor(plot_data_tissue$Biotype_binary, levels = c("protein-coding", "non-coding"))

stacked_bar_plot_tissue <- ggplot(plot_data_tissue, aes(x = tissue_labels, y = unique_genes, fill = Biotype_binary)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("protein-coding" = "#FFC20A", "non-coding" = "#0C7BDC")) +
  labs(title = '', x = "Tissue type", y = "Number of genes", fill = "Gene biotype") +
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits = c(0, 910)) +  # Set y-axis limits and labels
  theme_classic() +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, family = "Calibri", color = "black"),
    axis.text.y = element_text(size = 8, family = "Calibri", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 9, family = "Calibri", face = "bold", color = "black"),  # Enlarge plot title
    axis.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge axis labels
    legend.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge legend title
    legend.text = element_text(size = 8, family = "Calibri", color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(size = 0.3)) +
    theme(legend.key.height = unit(0.2, "mm"), legend.key.width = unit(2, "mm"))  # Adjust size of legend keys

stacked_bar_plot_tissue

ggsave(stacked_bar_plot_tissue, filename = "Figure_4/Figure_4A/Fig_4A_Barplot_pheno_associat_All_gene_count_by_tissue.pdf", width = 2.85, height = 2.18, device = cairo_pdf, dpi = 800)


# II.
# Figure 4B
# CAD novel candidate genes

# Select only novel CAD genes
novel_genes <- coloc_output[coloc_output$known_CAD_loci_gene_new_list == 'False', ]
length(unique(novel_genes$Gene.name))

# 4B - 1
# Gene freq by phenotype 
plot_data <- novel_genes
plot_data <- plot_data %>%
  dplyr::group_by(Trait_merged_groups, Biotype_binary) %>%
  dplyr::summarise(unique_genes = dplyr::n_distinct(Gene.name)) %>%
  dplyr::ungroup()

# Calculate the total number of unique genes per trait and sort by it
plot_data <- plot_data %>%
  group_by(Trait_merged_groups) %>%
  mutate(total_unique_genes = sum(unique_genes)) %>%
  ungroup() %>%
  arrange(desc(total_unique_genes))

trait_order <- plot_data %>%
  dplyr::group_by(Trait_merged_groups) %>%
  dplyr::summarise(total_unique_genes = sum(unique_genes)) %>%
  dplyr::arrange(desc(total_unique_genes)) %>%
  dplyr::pull(Trait_merged_groups)

# Ensure that the Trait_merged_groups levels are ordered from high to low based on the total number of unique genes
plot_data$Trait_merged_groups <- factor(plot_data$Trait_merged_groups, levels = trait_order)

# Create a barplot
plot_data$Biotype_binary <- factor(plot_data$Biotype_binary, levels = c("protein-coding", "non-coding"))

stacked_bar_plot <- ggplot(plot_data, aes(x = Trait_merged_groups, y = unique_genes, fill = Biotype_binary)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("protein-coding" = "#FFC20A", "non-coding" = "#0C7BDC")) +
  labs(title = '', x = "Phenotype", y = "Number of genes", fill = "Gene biotype") +
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 410)) +  # Set y-axis limits and labels
  theme_classic() +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, family = "Calibri", color = "black"),
    axis.text.y = element_text(size = 8, family = "Calibri", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 9, family = "Calibri", face = "bold", color = "black"),  # Enlarge plot title
    axis.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge axis labels
    legend.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge legend title
    legend.text = element_text(size = 8, family = "Calibri", color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(size = 0.3)) +
  
    theme(legend.key.height = unit(0.2, "mm"), legend.key.width = unit(2, "mm"))  # Adjust size of legend keys

stacked_bar_plot

ggsave(stacked_bar_plot, filename = "Figure_4/Figure_4A/Fig_4A_Barplot_pheno_associat_Novel_gene_count_by_pheno.pdf", width = 3.05, height = 2.88, device = cairo_pdf, dpi = 800)

# 4B - 2
# Gene count by tissue
plot_data_tissue <- novel_genes

# Group by tissue and count
plot_data_tissue <- plot_data_tissue %>%
  dplyr::group_by(tissue_labels, Biotype_binary) %>%
  dplyr::summarise(unique_genes = dplyr::n_distinct(Gene.name)) %>%
  dplyr::ungroup()

# Calculate the total number of unique genes per tissue and sort by it
plot_data_tissue <- plot_data_tissue %>%
  group_by(tissue_labels) %>%
  mutate(total_unique_genes = sum(unique_genes)) %>%
  ungroup() %>%
  arrange(desc(total_unique_genes))

tissue_order <- plot_data_tissue %>%
  dplyr::group_by(tissue_labels) %>%
  dplyr::summarise(total_unique_genes = sum(unique_genes)) %>%
  dplyr::arrange(desc(total_unique_genes)) %>%
  dplyr::pull(tissue_labels)

# Ensure that the tissue_labels levels are ordered from high to low based on the total number of unique genes
plot_data_tissue$tissue_labels <- factor(plot_data_tissue$tissue_labels, levels = tissue_order)

# Create a barplot
# Convert Biotype_binary to a factor with the desired order
plot_data_tissue$Biotype_binary <- factor(plot_data_tissue$Biotype_binary, levels = c("protein-coding", "non-coding"))

stacked_bar_plot_tissue <- ggplot(plot_data_tissue, aes(x = tissue_labels, y = unique_genes, fill = Biotype_binary)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("protein-coding" = "#FFC20A", "non-coding" = "#0C7BDC")) +
  labs(title = '', x = "Tissue type", y = "Number of genes", fill = "Gene biotype") +
  scale_y_continuous(breaks = seq(0, 400, by = 100), limits = c(0, 410)) +  # Set y-axis limits and labels
  theme_classic() +
    theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, family = "Calibri", color = "black"),
    axis.text.y = element_text(size = 8, family = "Calibri", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 9, family = "Calibri", face = "bold", color = "black"),  # Enlarge plot title
    axis.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge axis labels
    legend.title = element_text(size = 8, family = "Calibri", face = "bold", color = "black"),  # Enlarge legend title
    legend.text = element_text(size = 8, family = "Calibri", color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(size = 0.3)) +
    theme(legend.key.height = unit(0.2, "mm"), legend.key.width = unit(2, "mm"))  # Adjust size of legend keys

stacked_bar_plot_tissue

ggsave(stacked_bar_plot_tissue, filename = "Figure_4/Figure_4B/Fig_4B_Barplot_pheno_associat_Novel_gene_count_by_tissue.pdf", width = 2.85, height = 2.18, device = cairo_pdf, dpi = 800)
