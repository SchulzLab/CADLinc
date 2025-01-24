library(ggplot2)
library(ComplexHeatmap)
library(stylo) #for cosine distance
library(RColorBrewer)
library(circlize)
library(rcartocolor)

## there are 39 TFs with an odds-ratio higher 5 and more than 5 rSNVs which are also active (open) in their cell-type with high odds ratio
SIG_TFs = c(
"KLF12",
"ZIC4",
"KLF15",
"MBD2",
"ETV5::FOXO1",
"JUN",
"FOS::JUN",
"ZNF93",
"ZFX",
"MAZ",
"KLF7",
"KLF14",
"TFAP2C",
"ZIC5",
"NRF1",
"SP4",
"CTCFL",
"ZNF436",
"PATZ1",
"ZBTB33",
"FOSL2::JUN",
"TFEB",
"BATF3",
"FOS",
"CTCF",
"TFAP2E",
"FOS::JUND",
"KLF10",
"ZIC1",
"ZNF610",
"NFIX",
"FOSB",
"ZFP14",
"SP1",
"KLF16",
"TCFL5",
"THRB",
"ZNF707"
)

#input files
args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # hetamap odds ratio
output <- args[2] # output file gainLoss


## read in heatmap odds ratio 
data_ <- read.delim(inputFile, header = TRUE, sep = ",") # read in data for odds ratio

## fix issue with the dots in the header of a data frame
header = colnames(data_)
header = gsub("..", "::", header, fixed = TRUE )
colnames(data_) = header

print(head(data_))
celltype = as.matrix(data_[, 2, drop = FALSE]) #save row name

rownames(data_) = celltype #add row names
#remove info columns from data matrix
drops = c("celltype", "label")
data_plot = data_[, !(names(data_)) %in% drops]

print(ncol(data_plot))
data_plot = data_[, (names(data_)) %in% SIG_TFs] # only keep sig TFs
print(head(data_plot))
print(ncol(data_plot))

knownCAD= as.matrix(data_plot[1,])
newCAD= as.matrix(data_plot[2,])
targetGenes= as.matrix(data_plot[3,])
regSNVs= as.matrix(data_plot[4,])

## non-coding Genes, codingGenes and others (pseudo or TEC)
ncGenes = as.matrix(data_plot[5,])
cGenes = as.matrix(data_plot[6,])
otherGenes = as.matrix(data_plot[7,])

print("here")
data_plot = data_plot[-c(1,2,3,4,5,6,7),]
print(head(data_plot))

data_matrix = as.matrix(data_plot) #convert to matrix
header = colnames(data_matrix)
#print(length(header))

#plot heatmap
data_heatmap = t(data_plot)

#define clusting
hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist.cosine(x)

cluster_traits = hclustfunc(distfunc(data_heatmap))
ordering_traits = rownames(data_heatmap[cluster_traits$order,])

cluster_TFs = hclustfunc(distfunc(t(data_heatmap)))

helper = brewer.pal(9, "GnBu")
col_fun = colorRamp2(c(0, 5, 6, 7, 8, 9, 10, 15), c(helper[1], helper[3], helper[4] ,helper[5], helper[6], helper[7], helper[8], helper[9]))

#pdf(output) ## open output file
pdf(output, width=10) ## open output file
print(paste("minimal entry: ", min(data_heatmap), " maximal entry: ", max(data_heatmap), " number rows: ", nrow(data_heatmap), " number cols: ", ncol(data_heatmap), sep = ""))

row_bar = rowAnnotation(bar2 = anno_barplot(cbind(t(knownCAD), t(newCAD), t(targetGenes)), nc = 3, gp = gpar(fill = c("#E66100", "#5D3A9B", "#999999"), col ="grey"), bar_width = 0.85, border = FALSE, axis_param = list(gp=gpar(fontsize = 12),direction = "reverse" )), foo_empty = anno_empty(border = FALSE, width =unit(0.3, "cm")),  bar1 = anno_barplot(cbind(t(cGenes), t(ncGenes)), nc = 2, gp = gpar(fill = c("#FFC20A", "#0C7BDC"), col ="grey"), bar_width = 0.85, border = FALSE, axis_param = list(gp=gpar(fontsize = 12),direction = "reverse" )),  
		show_annotation_name = FALSE, 
		width = unit(2.5, "cm"), 
		show_legend = c(TRUE, FALSE, TRUE)
)

lgd = Legend(labels = c("known CAD GWAS genes","candidate CAD genes","others"), title = "Association to CAD", legend_gp = gpar(fill = c("#E66100", "#5D3A9B", "#999999"), fontsize = 8), title_gp = gpar(fontsize = 10), ncol = 3)
lgd2 = Legend(labels = c("protein-coding","non-coding RNA"), title = "Biotype", legend_gp = gpar(fill = c("#FFC20A", "#0C7BDC"), fontsize = 8), title_gp = gpar(fontsize = 10), ncol = 2)

ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm")

ht = Heatmap(data_heatmap, name = "?",
	heatmap_legend_param = list(title = "odds-ratio", title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8), at = c(0, 5,6,7, 8, 9, 10,  15), direction = "horizontal"),
	cluster_rows = as.dendrogram(cluster_traits),
	cluster_columns = as.dendrogram(cluster_TFs),
	col = col_fun, 
	show_column_names = TRUE,
	show_column_dend  = FALSE,
	column_title = "cell type", 
	column_title_gp = gpar(fontsize = 18),
	column_names_gp = gpar(fontsize = 10),
	column_title_side = "bottom", #"top"
	show_row_names = TRUE,
	show_row_dend = FALSE,
	row_title = "CAD-TF",
	row_title_gp = gpar(fontsize = 18),
	row_names_gp = gpar(fontsize = 9),
	row_title_side = "left", #"right", #"left"
	left_annotation = row_bar
)
draw(ht, annotation_legend_list = list(lgd,lgd2), merge_legend = TRUE, heatmap_legend_side = "top")
decorate_annotation("bar2", {grid.text("number of TF\ntarget genes", y = -0.17 , x = 1.2, just = "bottom", gp=gpar(fontsize = 12), rot = 0)}) ## change label of the annotation via decorate_annotation
dev.off()
