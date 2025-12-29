library(dplyr)
library(ggplot2)
library(plyr)
library(readxl)
library(tidyr)
library(tidyverse)
library(extrafont)
library(scales)
library(cowplot)
library(grid)
library(Cairo)
library(gridExtra)
library(patchwork)

coloc_output <- read.csv('results/gene_eQTL_GWAS_Joint_GeneBase/GWAS_extraction/Coloc/Final_merged_tables_GTEx_STARNET/Coloc_GTEx_STARNET_merged_eQTLs_HPP4_08_P_5e_8_Renamed_pheno_renamed_tissue.csv')

# 481 Unique CAD new genes that  have significant colocalization results
# 692 Unique CAD known genes that have significant colocalization results
# Count unique gene_ensembl values with 'False' in the 'known_CAD_loci_gene_new_list' column
false_count_known <- coloc_output %>%
  dplyr::filter(known_CAD_loci_gene_new_list == "False") %>%
  dplyr::summarise(unique_false = n_distinct(gene_ensembl))

# Count unique gene_ensembl values with 'True' in the 'known_CAD_loci_gene_new_list' column
true_count_new <- coloc_output %>%
  dplyr::filter(known_CAD_loci_gene_new_list == "True") %>%
  dplyr::summarise(unique_true = n_distinct(gene_ensembl))

# Print the results
false_count_known
true_count_new

# 1010 protein-coding and 163 non-coding genes
biotype_genes <- coloc_output[, c("Gene.name", "Biotype_binary")]
biotype_genes <- biotype_genes[!duplicated(biotype_genes), ]
table(biotype_genes$Biotype_binary)

# Select only novel CAD genes
novel_genes_coloc <- coloc_output[coloc_output$known_CAD_loci_gene_new_list == 'False', ]

# Filter dataframe based on 'Biotype_binary'
filtered_df_prot <- novel_genes_coloc[novel_genes_coloc$Biotype_binary == 'protein-coding', ]
filtered_df_lnc <- novel_genes_coloc[novel_genes_coloc$Biotype_binary == 'non-coding', ]

# 387 prot-coding novel CAD genes, 94 non-coding novel CAD genes
# 481 novel in total
length(unique(filtered_df_prot$Gene.name))
length(unique(filtered_df_lnc$Gene.name))




# I.
# Figure 4C
# Protein-coding

# Final selection
# There are way too many genes with 0 p-value as the lowest p-value
# So first, select all top genes with 0 p-value. Then count, how many traits they are associated with. 
# Step 1: Group by gene and trait, then select the row with the lowest p_value_GWAS in each group
  group_by(Gene.name, Trait) %>%
  slice_min(order_by = p_value_GWAS, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(p_value_GWAS)

# Step 2: Filter genes that have 0 as the lowest p-value
genes_with_zero_pvalue <- top_prot_cod_genes %>%
  dplyr::filter(p_value_GWAS == 0) %>%
  distinct(Gene.name)

# Step 3: Filter the original dataset to include only genes with zero p-value
filtered_genes_data <- filtered_df_prot %>%
  dplyr::filter(Gene.name %in% genes_with_zero_pvalue$Gene.name)

# Step 4: Count the number of different traits for each of those selected genes in the original data frame
trait_counts <- dplyr::filter(filtered_genes_data, Gene.name %in% genes_with_zero_pvalue$Gene.name) %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::summarize(distinct_traits_count = dplyr::n_distinct(Trait), .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(distinct_traits_count))

# Step 5: Select the top 10 genes with the highest number of traits
top_10_genes <- trait_counts %>%
  slice_head(n = 10)

selected_top_prot_coding <- filtered_df_prot %>%
  dplyr::filter(Gene.name %in% top_10_genes$Gene.name)

unique(selected_top_prot_coding$Gene.name)
write.csv(selected_top_prot_coding, 'Figure_4/Tables/Top_10_prot_coding_genes_with_all_significant_SNPs.csv', row.names = FALSE)

# Select top SNP for each gene in each combination with the lowest p-value in p_value_GWAS
top_snp_prot_cod <- selected_top_prot_coding %>%
  group_by(Gene.name, Trait_labels, tissue_labels) %>%
  slice(which.min(p_value_GWAS)) %>%
  ungroup()

top_snp_prot_cod$p_value_GWAS_log10 <- -log10(top_snp_prot_cod$p_value_GWAS)
write.csv(top_snp_prot_cod, 'Figure_4/Tables/Top_10_prot_coding_genes_with_TOP_SNPs.csv', row.names = FALSE)

# Change Inf log10 to the max value available in the df
top_snp_prot_cod <- top_snp_prot_cod %>%
  mutate(p_value_GWAS_log10 = ifelse(is.infinite(p_value_GWAS_log10), 323, p_value_GWAS_log10))

# Top tissue selection
# Compare the lowest p-value for each tissue. The one with the lowest - top tissue. If the values are the same, count the number of different traits for every tissue for the gene. The one with the largest number of different traits - top tissue
# Step 1: Calculate the minimum p-values for each gene and tissue
min_pvals_per_gene_tissue <- top_snp_prot_cod %>%
  dplyr::group_by(Gene.name, tissue_labels) %>%
  dplyr::summarise(
    min_pval = min(c(p_value_eQTL, p_value_GWAS), na.rm = TRUE),
    .groups = 'drop'
  )

# Step 2: Identify the tissue with the lowest p-value for each gene
lowest_pval_tissue <- min_pvals_per_gene_tissue %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::filter(min_pval == min(min_pval)) %>%
  dplyr::ungroup()

# Step 3: Count the number of unique 'Trait' values for each gene and tissue
unique_traits_count <- top_snp_prot_cod %>%
  dplyr::group_by(Gene.name, tissue_labels) %>%
  dplyr::summarise(unique_trait_count = n_distinct(Trait), .groups = 'drop')

# Step 4: Join the lowest p-value tissue information with the unique traits count
combined_df <- lowest_pval_tissue %>%
  dplyr::left_join(unique_traits_count, by = c("Gene.name", "tissue_labels"))
write.csv(combined_df, 'Figure_4/Tables/Top_10_prot_coding_genes_lowest_P_Trait_by_tissue_count.csv')

# Step 5: If there is a tie in p-values, select the tissue with the highest number of unique traits
top_tissue_per_gene <- combined_df %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::filter(unique_trait_count == max(unique_trait_count)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Step 6: Create the final data frame with the results
top_tissues <- top_tissue_per_gene %>%
  dplyr::select(Gene.name, tissue_labels)

top_tissues

# Filter the original dataframe to keep only the rows with the top tissues
top_snp_prot_cod_top_tissue <- top_snp_prot_cod %>%
  semi_join(top_tissues, by = c("Gene.name", "tissue_labels"))

# Change the Trait order
# Define the custom order for the levels
new_trait_order_prot_cod <- c("Nonalcoholic fatty liver disease", 'Type 2 diabetes', 'Two hour glucose test', 'HbA1c', 'Fasting insulin' , 'Fasting glucose', 
                     'White blood cell count', 'Red blood cell count', 'Platelet count', 'Neutrophil count', 'Monocyte count', 'Lymphocyte count', 'Eosinophil count', 'Basophil count',
                     'Heart failure', 'Resting heart rate', 'Peripheral arterial disease', 'Pulse pressure', 'Systolic blood pressure', 'Diastolic blood pressure', 
                     'Waist-hip ratio adj to BMI', 'Waist-hip ratio', 'Body mass index', 
                     'Triglycerides', 'Total cholesterol', 'Non-HDL cholesterol', 'HDL cholesterol', 'LDL cholesterol')

# Set the 'Trait_labels' column as a factor with the custom order
top_snp_prot_cod_top_tissue$Trait_labels <- factor(top_snp_prot_cod_top_tissue$Trait_labels, levels = new_trait_order_prot_cod)

# Add colors for top tissue for the gene
# Create a data frame with unique genes and their tissue labels
gene_colors <- top_snp_prot_cod_top_tissue %>%
  distinct(Gene.name, tissue_labels) %>%
  mutate(gene_color = ifelse(grepl('Artery', tissue_labels), 'darkred', 
                             ifelse(grepl('Adipose tissue', tissue_labels), 'darkorange', 
                                    ifelse(grepl('Blood', tissue_labels), '#F8a99d',
                                          ifelse(grepl('Heart', tissue_labels), '#B378D3', 
                                           ifelse(grepl('Liver', tissue_labels), '#66c2a5',
                                                  ifelse(grepl('Skeletal muscle', tissue_labels), '#3288bd',
                                                         ifelse(grepl('Spleen', tissue_labels), '#1d2951',
                                                                'black'))))))))

# Convert gene_colors to a named vector for easy lookup
gene_colors_vec <- setNames(gene_colors$gene_color, gene_colors$Gene.name)
gene_colors_vec

# Add the conservation part
# Convert to 'yes'/'no'
top_snp_prot_cod_top_tissue$Conserved_in_mouse <- ifelse(top_snp_prot_cod_top_tissue$Conserved.in.mouse, "yes", "no")
top_snp_prot_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'
top_snp_prot_cod_top_tissue$Conserved_in_mouse <- factor(top_snp_prot_cod_top_tissue$Conserved_in_mouse,
                                                         levels = c("no", "yes"))

top_snp_prot_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'
head(top_snp_prot_cod_top_tissue)

# Protein-coding
# Combined pdf figure
heatmap_plot_prot_cod <- ggplot(top_snp_prot_cod_top_tissue, aes(x = Gene.name, y = Trait_labels)) +
  geom_tile(aes(fill = p_value_GWAS_log10), size = 0.2) + 
  scale_fill_gradientn(name = "-log10(p-value)",
                       colours = c("white", 'blue', "#00008B", "#010048"), 
                       values = scales::rescale(c(0, 100, 200, 323), to = c(0, 1)),
                       breaks = c(0, 100, 200, 300), 
                       labels = c("0", "100", '200', '300'),
                       na.value = 'white',
                       limits = c(0, 323)) +
  theme_classic() +
  labs(title = "", x = "", y = "Phenotype") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),  
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),  
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),  
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 7, face = "bold", color = "black"),
    axis.line = element_line(size = 0), 
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),  
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0)  # Remove extra margin
  ) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_prot_cod_top_tissue$Trait_labels)) + 0.5), color = "grey", size = 0.1) + 
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_prot_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1)

# Second heatmap
heatmap_conservation_prot_cod <- ggplot(top_snp_prot_cod_top_tissue, aes(x = Gene.name, y = Supportive_evidence)) +
  geom_tile(aes(fill = Conserved_in_mouse), size = 0.2) +
  scale_fill_manual(name = "Supportive evidence",
                    values = c("yes" =  "#ffcc66",  "no" = "white"),
                    labels = c("no", "yes"),
                    drop = FALSE) +

  theme_classic() +
  labs(title = "", x = "", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold.italic",
                               color = gene_colors_vec[levels(factor(top_snp_prot_cod_top_tissue$Gene.name))]),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 7, color = "black", angle = 0),
    axis.line = element_line(size = 0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),  # Remove extra margin
    aspect.ratio = length(unique(top_snp_prot_cod_top_tissue$Supportive_evidence)) / length(unique(top_snp_prot_cod_top_tissue$Gene.name))
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_prot_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_prot_cod_top_tissue$Conserved_in_mouse))), color = "grey", size = 0.1)

# Extract legends
legend1 <- get_legend(heatmap_plot_prot_cod)
legend2 <- get_legend(heatmap_conservation_prot_cod)

# Remove legends from original plots
heatmap_plot_prot_cod <- heatmap_plot_prot_cod + theme(legend.position = "none")
heatmap_conservation_prot_cod <- heatmap_conservation_prot_cod + theme(legend.position = "none")

# Create a dummy plot for the additional tissue type legend
dummy_data <- data.frame(
  category = c("Artery", "Adipose", "Blood", 'Heart', 'Liver', 'SKLM', 'Spleen'),
  value = c(1, 1, 1, 1, 1, 1, 1)
)

tissue_legend_plot <- ggplot(dummy_data, aes(x = category, fill = category)) +
  geom_bar() +
  scale_fill_manual(name = "Tissue type",
                    values = c("Artery" = 'darkorange', "Adipose" = "darkred", "Blood" = '#F8a99d', "Heart" = '#B378D3', "Liver" = '#66c2a5', "SKLM" = '#3288bd', 'Spleen' = '#1d2951'),
                    labels = c("Artery", "Adipose", 'Blood', 'Heart', 'Liver', 'SKLM', 'Spleen')) +
  theme_void() +
  theme(
    legend.title = element_text(size = 7, color = "black", face = "bold"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black")
  ) 

# Extract the legend
legend3 <- get_legend(tissue_legend_plot)

combined_legends <- plot_grid(
  legend1, legend2, legend3, 
  ncol = 1, 
  align = 'v',
  rel_heights = c(1, 1, 1) 
)

# Align heatmaps
aligned_plots <- cowplot::align_plots(heatmap_plot_prot_cod, heatmap_conservation_prot_cod, align = 'v', axis = 'b')

# Combine the aligned heatmaps
combined_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]], 
  ncol = 1, 
  rel_heights = c(1, 1),
  align = 'v'
)

combined_plot

# Combine the plots and legends
final_plot <- plot_grid(
  combined_plot, 
  combined_legends, 
  ncol = 2, 
  rel_widths = c(1, 1)
)

final_plot

# Save the final plot
file_path <- "Figure_4/Figure_4C/Prot_coding_combined_Figure.pdf"
CairoPDF(file_path, width = 43, height = 53, dpi = 600)
print(final_plot)
dev.off()




# II.
# Figure 4D
# Non-coding

# Final selection
# Group by gene and trait, then select the row with the lowest p_value_GWAS in each group
top_non_cod_genes <- filtered_df_lnc %>%
  dplyr::group_by(gene_ensembl, Trait) %>%
  dplyr::filter(p_value_GWAS == min(p_value_GWAS)) %>%
  ungroup() %>%
  arrange(p_value_GWAS)

top_10_genes <- top_non_cod_genes %>%
  distinct(gene_ensembl, .keep_all = TRUE) %>%
  slice_head(n = 10)

selected_top_non_coding <- filtered_df_lnc %>%
  dplyr::filter(Gene.name %in% top_10_genes$Gene.name)

unique(selected_top_non_coding$Gene.name)
write.csv(selected_top_non_coding, 'Figure_4/Tables/Top_10_non_coding_genes_with_all_significant_SNPs.csv', row.names = FALSE)

# Select top SNP for each gene in each combination with the lowest p-value in p_value_GWAS
top_snp_non_cod <- selected_top_non_coding %>%
  group_by(Gene.name, Trait_labels, tissue_labels) %>%
  slice(which.min(p_value_GWAS)) %>%
  ungroup()

top_snp_non_cod$p_value_GWAS_log10 <- -log10(top_snp_non_cod$p_value_GWAS)
write.csv(top_snp_non_cod, 'Figure_4/Tables/Top_10_non_coding_genes_with_TOP_SNPs.csv', row.names = FALSE)

top_snp_non_cod <- top_snp_non_cod %>%
  mutate(p_value_GWAS_log10 = ifelse(is.infinite(p_value_GWAS_log10), 316, p_value_GWAS_log10))

# Top tissue selection
# Step 1: Calculate the minimum p-values for each gene and tissue
min_pvals_per_gene_tissue <- top_snp_non_cod %>%
  dplyr::group_by(Gene.name, tissue_labels) %>%
  dplyr::summarise(
    min_pval = min(c(p_value_eQTL, p_value_GWAS), na.rm = TRUE),
    .groups = 'drop'
  )

# Step 2: Identify the tissue with the lowest p-value for each gene
lowest_pval_tissue <- min_pvals_per_gene_tissue %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::filter(min_pval == min(min_pval)) %>%
  dplyr::ungroup()

# Step 3: Count the number of unique 'Trait' values for each gene and tissue
unique_traits_count <- top_snp_non_cod %>%
  dplyr::group_by(Gene.name, tissue_labels) %>%
  dplyr::summarise(unique_trait_count = n_distinct(Trait), .groups = 'drop')

# Step 4: Join the lowest p-value tissue information with the unique traits count
combined_df <- lowest_pval_tissue %>%
  dplyr::left_join(unique_traits_count, by = c("Gene.name", "tissue_labels"))

write.csv(combined_df, 'Figure_4/Tables/Top_10_non_coding_genes_lowest_P_Trait_by_tissue_count.csv')

# Step 5: If there is a tie in p-values, select the tissue with the highest number of unique traits
top_tissue_per_gene <- combined_df %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::filter(unique_trait_count == max(unique_trait_count)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Step 6: Create the final data frame with the results
top_tissues <- top_tissue_per_gene %>%
  dplyr::select(Gene.name, tissue_labels)

top_tissues

# Filter the original dataframe to keep only the rows with the top tissues
top_snp_non_cod_top_tissue <- top_snp_non_cod %>%
  semi_join(top_tissues, by = c("Gene.name", "tissue_labels"))

# Change the Trait order
# Define the custom order for the levels
new_trait_order_prot_cod <- c("Nonalcoholic fatty liver disease", 'Type 2 diabetes', 'Two hour glucose test', 'HbA1c', 'Fasting insulin' , 'Fasting glucose', 
                     'White blood cell count', 'Red blood cell count', 'Platelet count', 'Neutrophil count', 'Monocyte count', 'Lymphocyte count', 'Eosinophil count', 'Basophil count',
                     'Heart failure', 'Resting heart rate', 'Peripheral arterial disease', 'Pulse pressure', 'Systolic blood pressure', 'Diastolic blood pressure', 
                     'Waist-hip ratio adj to BMI', 'Waist-hip ratio', 'Body mass index', 
                     'Triglycerides', 'Total cholesterol', 'Non-HDL cholesterol', 'HDL cholesterol', 'LDL cholesterol')

# Set the 'Trait_labels' column as a factor with the custom order
top_snp_non_cod_top_tissue$Trait_labels <- factor(top_snp_non_cod_top_tissue$Trait_labels, levels = new_trait_order_prot_cod)

# Add colors for top tissue for the gene
# Create a data frame with unique genes and their tissue labels
gene_colors_non_cod <- top_snp_non_cod_top_tissue %>%
  distinct(Gene.name, tissue_labels) %>%
  mutate(gene_color = ifelse(grepl('Artery', tissue_labels), 'darkred', 
                             ifelse(grepl('Adipose tissue', tissue_labels), 'darkorange', 
                                    ifelse(grepl('Blood', tissue_labels), '#F8a99d',
                                          ifelse(grepl('Heart', tissue_labels), '#B378D3', 
                                           ifelse(grepl('Liver', tissue_labels), '#66c2a5',
                                                  ifelse(grepl('Skeletal muscle', tissue_labels), '#3288bd',
                                                         ifelse(grepl('Spleen', tissue_labels), '#1d2951',
                                                                'black'))))))))

# Convert gene_colors_non_cod to a named vector for easy lookup
gene_colors_non_cod_vec <- setNames(gene_colors_non_cod$gene_color, gene_colors_non_cod$Gene.name)
gene_colors_non_cod_vec

# Add the conservation part
## Convert to 'yes'/'no'
top_snp_non_cod_top_tissue$Conserved_in_mouse <- ifelse(top_snp_non_cod_top_tissue$Conserved.in.mouse, "yes", "no")
top_snp_non_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'

top_snp_non_cod_top_tissue$Conserved_in_mouse <- factor(top_snp_non_cod_top_tissue$Conserved_in_mouse,
                                                         levels = c("no", "yes"))
top_snp_non_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'
head(top_snp_non_cod_top_tissue)

# Non-coding
# Combined pdf figure
heatmap_plot_non_cod <- ggplot(top_snp_non_cod_top_tissue, aes(x = Gene.name, y = Trait_labels)) +
  geom_tile(aes(fill = p_value_GWAS_log10), size = 0.2) + 
  scale_fill_gradientn(name = "-log10(p-value)",
                       colours = c("white", 'blue', "#00008B", "#010048"), 
                       values = scales::rescale(c(0, 100, 200, 320), to = c(0, 1)),
                       breaks = c(0, 100, 200, 300), 
                       labels = c("0", "100", '200', '300'),
                       na.value = 'white',
                       limits = c(0, 316)) +
  theme_classic() +
  labs(title = "", x = "", y = "Phenotype") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),  
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),  
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),  
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 7, face = "bold", color = "black"),
    axis.line = element_line(size = 0), 
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),  
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0)  # Remove extra margin
  ) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_non_cod_top_tissue$Trait_labels)) + 0.5), color = "grey", size = 0.1) + 
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_non_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1)

# Second heatmap
heatmap_conservation_non_cod <- ggplot(top_snp_non_cod_top_tissue, aes(x = Gene.name, y = Supportive_evidence)) +
  geom_tile(aes(fill = Conserved_in_mouse), size = 0.2) +
  scale_fill_manual(name = "Supportive evidence",
                    values = c("yes" = "#ffcc66", "no" = "white"),
                    labels = c("no", "yes"),
                    drop = FALSE) +
  theme_classic() +
  labs(title = "", x = "", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold.italic",
                               color = gene_colors_non_cod_vec[levels(factor(top_snp_non_cod_top_tissue$Gene.name))]),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 7, color = "black", angle = 0),
    axis.line = element_line(size = 0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),  # Remove extra margin
    aspect.ratio = length(unique(top_snp_non_cod_top_tissue$Supportive_evidence)) / length(unique(top_snp_non_cod_top_tissue$Gene.name))
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_non_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_non_cod_top_tissue$Conserved_in_mouse)) + 0.5), color = "grey", size = 0.1)

# Extract legends
legend1 <- get_legend(heatmap_plot_non_cod)
legend2 <- get_legend(heatmap_conservation_non_cod)

# Remove legends from original plots
heatmap_plot_non_cod <- heatmap_plot_non_cod + theme(legend.position = "none")
heatmap_conservation_non_cod <- heatmap_conservation_non_cod + theme(legend.position = "none")

# Create a dummy plot for the additional tissue type legend
dummy_data <- data.frame(
  category = c("Artery", "Adipose", "Blood", 'Heart', 'Liver', 'SKLM', 'Spleen'),
  value = c(1, 1, 1, 1, 1, 1, 1)
)

tissue_legend_plot <- ggplot(dummy_data, aes(x = category, fill = category)) +
  geom_bar() +
  scale_fill_manual(name = "Tissue type",
                    values = c("Artery" = 'darkorange', "Adipose" = "darkred", "Blood" = '#F8a99d', "Heart" = '#B378D3', "Liver" = '#66c2a5', "SKLM" = '#3288bd', 'Spleen' = '#1d2951'),
                    labels = c("Artery", "Adipose", 'Blood', 'Heart', 'Liver', 'SKLM', 'Spleen')) +
  theme_void() +
  theme(
    legend.title = element_text(size = 7, color = "black", face = "bold"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black")
  ) 

# Extract the legend
legend3 <- get_legend(tissue_legend_plot)

combined_legends <- plot_grid(
  legend1, legend2, legend3, 
  ncol = 1, 
  align = 'v',
  rel_heights = c(1, 1, 1)  # Adjust heights to bring legends closer
)

# Align heatmaps
aligned_plots <- cowplot::align_plots(heatmap_plot_non_cod, heatmap_conservation_non_cod, align = 'v', axis = 'b')

# Combine the aligned heatmaps
combined_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]], 
  ncol = 1, 
  rel_heights = c(1, 1),  
  align = 'v'
)

combined_plot

# Combine the plots and legends
final_plot <- plot_grid(
  combined_plot, 
  combined_legends, 
  ncol = 2, 
  rel_widths = c(1, 1)
)

final_plot

# Save the final plot
file_path <- "Figure_4/Figure_4D/Non_coding_combined_Figure.pdf"
CairoPDF(file_path, width = 43, height = 55, dpi = 600)
print(final_plot)
dev.off()




# III.
# Figure 4E
# Novel non-coding CAD genes, only conserved in mouse 

# Only 6 conserved in mouse lncRNAs have significant results
conserved_gene_count <- filtered_df_lnc %>%
  dplyr::filter(Conserved.in.mouse == 'True') %>%
  dplyr::summarize(count = dplyr::n_distinct(Gene.name))

conserved_gene_count

conserved_lnc <- filtered_df_lnc %>%
  dplyr::filter(Conserved.in.mouse == 'True')

conserved_lnc <- conserved_lnc %>%
  dplyr::filter(known_CAD_loci_gene_new_list == 'False')

# Count the number of different traits for each of those selected genes
trait_counts <- conserved_lnc %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::summarize(distinct_traits_count = dplyr::n_distinct(Trait), .groups = 'drop') %>%
  dplyr::arrange(dplyr::desc(distinct_traits_count))

write.csv(conserved_lnc, 'Figure_4/Tables/Novel_conserved_non_coding_genes_with_all_significant_SNPs.csv', row.names = FALSE)
conserved_lnc

# Select top SNP for each gene in each combination with the lowest p-value in p_value_GWAS
top_snp_non_cod_conserved <- conserved_lnc %>%
  group_by(Gene.name, Trait_labels, tissue_labels) %>%
  slice(which.min(p_value_GWAS)) %>%
  ungroup()

top_snp_non_cod_conserved$p_value_GWAS_log10 <- -log10(top_snp_non_cod_conserved$p_value_GWAS)
write.csv(top_snp_non_cod_conserved, 'Figure_4/Tables/Novel_conserved_non_coding_genes_with_TOP_SNPs.csv', row.names = FALSE)

# Top tissue selection
# Step 1: Count the number of unique 'Trait' values for each gene and tissue
unique_traits_count <- top_snp_non_cod_conserved %>%
  group_by(Gene.name, tissue_labels) %>%
  dplyr::summarise(unique_trait_count = n_distinct(Trait), .groups = "drop")

# Step 2: For each gene, keep the tissue with the highest number of unique traits
combined_df <- unique_traits_count %>%
  group_by(Gene.name) %>%
  dplyr::filter(unique_trait_count == max(unique_trait_count)) %>%
  slice(1) %>%   
  ungroup() %>%
  dplyr::select(Gene.name, tissue_labels, unique_trait_count)

write.csv(combined_df, 'Figure_4/Tables/Novel_conserved_non_coding_genes_lowest_P_Trait_by_tissue_count.csv')

# Step 3: Create the final data frame with the results
top_tissues_conserved <- combined_df %>%
  dplyr::select(Gene.name, tissue_labels)
top_tissues_conserved

# Filter the original dataframe to keep only the rows with the top tissues
top_snp_conserved_non_cod_top_tissue <- top_snp_non_cod_conserved %>%
  semi_join(top_tissues_conserved, by = c("Gene.name", "tissue_labels"))

top_snp_conserved_non_cod_top_tissue <- top_snp_conserved_non_cod_top_tissue %>%
  mutate(p_value_GWAS_log10 = ifelse(is.infinite(p_value_GWAS_log10), 157, p_value_GWAS_log10))

# Change the Trait order
# Define the custom order for the levels
new_trait_order <- c('Type 2 diabetes', 'Two hour glucose test', 'HbA1c', 'Fasting insulin' , 'Fasting glucose', 
                     'White blood cell count', 'Red blood cell count', 'Platelet count', 'Neutrophil count', 'Monocyte count', 'Lymphocyte count', 'Eosinophil count', 'Basophil count',
                     'Heart failure', 'Resting heart rate', 'Peripheral arterial disease', 'Pulse pressure', 'Systolic blood pressure', 'Diastolic blood pressure', 
                     'Waist-hip ratio adj to BMI', 'Waist-hip ratio', 'Body mass index', 
                     'Triglycerides', 'Total cholesterol', 'Non-HDL cholesterol', 'HDL cholesterol', 'LDL cholesterol')

# Set the 'Trait_labels' column as a factor with the custom order
top_snp_conserved_non_cod_top_tissue$Trait_labels <- factor(top_snp_conserved_non_cod_top_tissue$Trait_labels, levels = new_trait_order)

# Add colors for top tissue for the gene
# Create a data frame with unique genes and their tissue labels
gene_colors_conserved <- top_snp_conserved_non_cod_top_tissue %>%
  distinct(Gene.name, tissue_labels) %>%
  mutate(gene_color = ifelse(grepl('Artery', tissue_labels), 'darkred', 
                             ifelse(grepl('Adipose tissue', tissue_labels), 'darkorange', 
                                    ifelse(grepl('Blood', tissue_labels), '#F8a99d',
                                          ifelse(grepl('Heart', tissue_labels), '#B378D3', 
                                           ifelse(grepl('Liver', tissue_labels), '#66c2a5',
                                                  ifelse(grepl('Skeletal muscle', tissue_labels), '#3288bd',
                                                         ifelse(grepl('Spleen', tissue_labels), '#1d2951',
                                                                'black'))))))))

# Convert gene_colors_non_cod to a named vector for easy lookup
gene_colors_conserved_vec <- setNames(gene_colors_conserved$gene_color, gene_colors_conserved$Gene.name)

gene_colors_conserved_vec

# Add the conservation prt
# Convert to 'yes'/'no'
top_snp_conserved_non_cod_top_tissue$Conserved_in_mouse <- ifelse(top_snp_conserved_non_cod_top_tissue$Conserved.in.mouse, "yes", "no")
top_snp_conserved_non_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'

top_snp_conserved_non_cod_top_tissue$Conserved_in_mouse <- factor(top_snp_conserved_non_cod_top_tissue$Conserved_in_mouse,
                                                         levels = c("no", "yes"))

top_snp_conserved_non_cod_top_tissue$Supportive_evidence <- 'Conserved in mouse'
head(top_snp_conserved_non_cod_top_tissue)

# Non-coding, novel and conserved
# Combined pdf figure
heatmap_plot_non_cod <- ggplot(top_snp_conserved_non_cod_top_tissue, aes(x = Gene.name, y = Trait_labels)) +
  geom_tile(aes(fill = p_value_GWAS_log10), size = 0.2) + 
  scale_fill_gradientn(name = "-log10(p-value)",
                       colours = c("white", 'blue', "#00008B", "#010048"), 
                       values = scales::rescale(c(0, 100, 200, 320), to = c(0, 1)),
                       breaks = c(0, 100, 200, 300), 
                       labels = c("0", "100", '200', '300'),
                       na.value = 'white',
                       limits = c(0, 316)) +
  theme_classic() +
  labs(title = "", x = "", y = "Phenotype") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),  
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),  
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),  
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 7, face = "bold", color = "black"),
    axis.line = element_line(size = 0), # was 0.2-0.3
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),  
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0)  # Remove extra margin
  ) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_conserved_non_cod_top_tissue$Trait_labels)) + 0.5), color = "grey", size = 0.1) + 
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_conserved_non_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1)

# Second heatmap
heatmap_conservation_non_cod <- ggplot(top_snp_conserved_non_cod_top_tissue, aes(x = Gene.name, y = Supportive_evidence)) +
  geom_tile(aes(fill = Conserved_in_mouse), size = 0.2) +
  scale_fill_manual(name = "Supportive evidence",
                    values = c("yes" = "#ffcc66", "no" = "white"),
                    labels = c("no", "yes"),
                    drop = FALSE) +
  theme_classic() +
  labs(title = "", x = "", y = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9, face = "bold", color = "black"),
    legend.title = element_text(size = 7, color = "black", face = "bold", hjust = 0),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold.italic",
                               color = gene_colors_conserved_vec[levels(factor(top_snp_conserved_non_cod_top_tissue$Gene.name))]),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 7, color = "black", angle = 0),
    axis.line = element_line(size = 0),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0),  # Remove extra margin
    aspect.ratio = length(unique(top_snp_conserved_non_cod_top_tissue$Supportive_evidence)) / length(unique(top_snp_conserved_non_cod_top_tissue$Gene.name))
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(top_snp_conserved_non_cod_top_tissue$Gene.name)) + 0.5), color = "grey", size = 0.1) +
  geom_hline(yintercept = seq(0.5, length(unique(top_snp_conserved_non_cod_top_tissue$Conserved_in_mouse))), color = "grey", size = 0.1)

# Extract legends
legend1 <- get_legend(heatmap_plot_non_cod)
legend2 <- get_legend(heatmap_conservation_non_cod)

# Remove legends from original plots
heatmap_plot_non_cod <- heatmap_plot_non_cod + theme(legend.position = "none")
heatmap_conservation_non_cod <- heatmap_conservation_non_cod + theme(legend.position = "none")

# Create a dummy plot for the additional tissue type legend
dummy_data <- data.frame(
  category = c("Artery", "Adipose", "Blood", 'Heart', 'Liver', 'SKLM', 'Spleen'),
  value = c(1, 1, 1, 1, 1, 1, 1)
)

tissue_legend_plot <- ggplot(dummy_data, aes(x = category, fill = category)) +
  geom_bar() +
  scale_fill_manual(name = "Tissue type",
                    values = c("Artery" = 'darkorange', "Adipose" = "darkred", "Blood" = '#F8a99d', "Heart" = '#B378D3', "Liver" = '#66c2a5', "SKLM" = '#3288bd', 'Spleen' = '#1d2951'),
                    labels = c("Artery", "Adipose", 'Blood', 'Heart', 'Liver', 'SKLM', 'Spleen')) +
  theme_void() +
  theme(
    legend.title = element_text(size = 7, color = "black", face = "bold"),
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.3, "cm"),
    legend.key = element_rect(color = "black")
  ) 

# Extract the legend
legend3 <- get_legend(tissue_legend_plot)

combined_legends <- plot_grid(
  legend1, legend2, legend3, 
  ncol = 1, 
  align = 'v',
  rel_heights = c(1, 1, 1)  # Adjust heights to bring legends closer
)

# Align heatmaps
aligned_plots <- cowplot::align_plots(heatmap_plot_non_cod, heatmap_conservation_non_cod, align = 'v', axis = 'b')

# Combine the aligned heatmaps
combined_plot <- plot_grid(
  aligned_plots[[1]], 
  aligned_plots[[2]], 
  ncol = 1, 
  rel_heights = c(1, 1),  # Adjust heights as necessary
  align = 'v'
)

combined_plot

# Combine the plots and legends
final_plot <- plot_grid(
  combined_plot, 
  combined_legends, 
  ncol = 2, 
  rel_widths = c(1, 1) 
)

final_plot

# Save the final plot
file_path <- "Figure_4/Figure_4E/Novel_conserved_non_coding_combined_Figure.pdf"
CairoPDF(file_path, width = 33, height = 47, dpi = 600)
print(final_plot)
dev.off()
