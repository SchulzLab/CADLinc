library(ggplot2)
library(cowplot)
library(ggrepel)

#input files
args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # TFs with high eQTL evidence  ../figure_tables/TFs_high_eQTL_evidence.txt
inputFile2 <-args[2] # functional enrichment for the TFs from the previous file ../figure_tables/TFs_functionalEnrichment_detailedPathways_average.txt 
outputFile<- args[3] 

data <- read.delim(inputFile, header = TRUE, sep = "\t") # read in data for odds ratio
print(head(data))

## Dot plot for eQTL fraction
p1 = ggplot(data, aes(x=TF, y = fraction)) + 
	geom_point(aes(size = rSNVs)) +
	scale_x_discrete(limits = c("NFIX","THRB", "CTCF", "PATZ1","TCFL5","KLF16", "TFEB", "ZNF610", "TFAP2E", "CTCFL")) + 
	coord_flip() + 
	xlab("CAD-TFs") +
	ylab("") + 
	ylim(0, 1.0) + 
	theme_classic() +
	theme(panel.grid.major.y = element_line(color = "gray", linewidth = 0.2),  panel.background = element_rect(fill = "white", colour = "grey"), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 12), legend.position = 'none')#+


##TODO: add plot of functional enrichment and combine them 
data2 <- read.delim(inputFile2, header = TRUE, sep = "\t") # read in data for odds ratio

p2 = ggplot(data2, aes(x=enrichedPathway, y = TF, colour = averagePvalue)) + 
	geom_point( aes( shape = source, size = averageFraction)) +
	scale_y_discrete(limits = c("NFIX","THRB", "CTCF", "PATZ1","TCFL5","KLF16", "TFEB", "ZNF610", "TFAP2E", "CTCFL")) + 
	scale_shape_manual("source", values = c(15, 16, 17, 18)) + 
	scale_size_continuous("fraction", range = c(4, 8)) + 
	scale_color_continuous("p-value") + 
	scale_x_discrete(labels = c(
		"Immune response and inflammation" = "Immune response\nand inflammation",
		"Lipid metabolism" = "Lipid metabolism", 
		"Metabolite transportation" = "Metabolite transportation", 
		"Innervation and synaptic transmission" = "Innervation and\nsynaptic transmission",
		"TGF beta signaling pathway" = "TGF beta signaling pathway", 
		"Vascular remodelling" = "Vascular remodelling"
)) + 
	xlab("enriched pathways") +
	#ylab("CAD-TFs") + 
	theme_classic() +
	theme(panel.grid.major.y = element_line(color = "gray", linewidth = 0.2),  panel.background = element_rect(fill = "white", colour = "grey"), axis.title.x = element_text(size = 14), axis.title.y = element_blank(), axis.text.x = element_text(size = 12, angle = -45, vjust = 0.8, hjust=0), axis.text.y = element_blank(), legend.key.height = unit(0.5, "cm"))#+

plot_grid(p1, p2, ncol = 2, align = 'h', rel_widths = c(1, 2))
ggsave(outputFile)




