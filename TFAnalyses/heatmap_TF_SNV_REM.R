library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#input files
args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # hetamap odds ratio
output <- args[2] # output file gainLoss

## read in heatmap 
data_ <- read.delim(inputFile, header = TRUE, sep = "\t") # read in data
print(head(data_))

## change dataframe colnames to get avoid the dots in the header
header = colnames(data_) ## all specical characters are set to dot
header = gsub(".", " ", header, fixed = TRUE) # replace . by space
header = gsub("TF SNV", "TF-SNV", header, fixed = TRUE) # add - 
header = gsub("CD4  T cell","CD4+ T cell", header, fixed = TRUE) ## add + 
header = gsub("CD8  T cell","CD8+ T cell", header, fixed = TRUE) ## add + 
#set new colmanes 
colnames(data_) = header

# get data for annotation bars for heatmap
geneNames = as.matrix(data_[, "Gene name", drop = FALSE]) #save row name
mergedRems = as.matrix(data_[, "Number merged CREs with  1  FDR TF-SNV from cell types with open promoter", drop = FALSE]) #save number of merged rem across celltypes
expression = as.matrix(data_[, "Expression breadth", drop = FALSE]) #save expression values
biotype = as.matrix(data_[, "Biotype", drop = FALSE])  ## collect the biotype which is either protein_coding or lncRNA
biotype = gsub("protein_coding", "protein-coding", biotype, fixed = TRUE)
biotype = gsub("lncRNA", "non-coding RNA", biotype, fixed = TRUE)
biotype = gsub("miRNA", "non-coding RNA", biotype, fixed = TRUE)

evidence = as.matrix(data_[, "known CAD GWAS gene", drop = FALSE]) #save info for evidence
## replace the strings
evidence = gsub("True", "yes", evidence, fixed = TRUE)
evidence = gsub("False", "no", evidence, fixed = TRUE)
coexpressed = as.matrix(data_[, "Co expressed with CAD GWAS genes", drop = FALSE]) # save whether a gene is co-expressed or not with a CAD loci
print(coexpressed)
coexpressed = gsub("Co-expressed with CAD GWAS genes", "yes", coexpressed, fixed = TRUE)
coexpressed = gsub("True", "yes", coexpressed, fixed = TRUE)
coexpressed = gsub("False", "no", coexpressed, fixed = TRUE)
coexpressed = gsub("insufficient co-expressed protein coding genes", "no", coexpressed, fixed = TRUE)
coexpressed = gsub("insufficient expression data", "no", coexpressed, fixed = TRUE)

coloc_gtex = as.matrix(data_[, "Coloc GTEx", drop = FALSE])
coloc_gtex = gsub("True", "yes", coloc_gtex, fixed = TRUE)
coloc_gtex = gsub("False", "no", coloc_gtex, fixed = TRUE)

coloc_starnet = as.matrix(data_[, "Coloc STARNET", drop = FALSE])
coloc_starnet = gsub("True", "yes", coloc_starnet, fixed = TRUE)
coloc_starnet = gsub("False", "no", coloc_starnet, fixed = TRUE)

conserved = as.matrix(data_[, "Conserved in mouse", drop = FALSE])
print(conserved)
conserved = gsub("True", "yes", conserved, fixed = TRUE)
conserved = gsub("False", "no", conserved, fixed = TRUE)

ensemblID = as.matrix(data_[, "Ensembl ID", drop = FALSE])
rownames(data_) = geneNames #add row names

#remove info columns from data matrix
drops = c("Ensembl ID","Gene name", "known CAD GWAS gene", "Coloc GTEx", "Coloc STARNET", "found via", "Schnitzler  CAD_V2G2P", "Schnitzler  PoPS", "Schnitzler  Hodonsky_eQTL_colocalization", "Schnitzler  Hodonsky_sQTL_colocalization", "Schnitzler  OpenTargetL2G","Schnitzler  van_der_Harst_2018", "Co expressed with CAD GWAS genes", "Expression breadth","Biotype","General biotype", "Number  1  FDR TF-SNVs" , "Chr","Start", "End", "Strand", "Cell types with open promoter", "Conserved in mouse", "Cell types with open promoter and  2 CREs with non LD TF-SNVs", "Number  1  FDR TFs SNVs", "Number affected TFs", "Number merged CREs with  1  FDR TF-SNV from cell types with open promoter", "GATES FDR")


print(ncol(data_))
data_plot = data_[, !(names(data_)) %in% drops]
print(ncol(data_plot))
print(head(data_plot))

data_heatmap = as.matrix(data_plot) #convert to matrix
header = colnames(data_heatmap)
print(length(header))
print(head(data_heatmap))

#define clusting
hclustfunc <- function(x) hclust(x, method="ward.D")
distfunc <- function(x) dist.cosine(x)

## color definition
helper = brewer.pal(9, "GnBu")
blues = brewer.pal(9, "Blues")
col_fun = colorRamp2(c(0, 1, 2, 3, 4, 5, 6,7, 8), c(helper[1], helper[3], helper[4] ,helper[5], helper[6], helper[7], helper[8], helper[9], blues[9]))

pdf(output, width=18) ## open output file

print(paste("minimal entry: ", min(data_heatmap), " maximal entry: ", max(data_heatmap), " number rows: ", nrow(data_heatmap), " number cols: ", ncol(data_heatmap), sep = ""))

## define colorcode for row_annotations 

oranges = brewer.pal(9, "Oranges")
col_fun_expression = colorRamp2(c(0,0.5, 1), c(oranges[1], oranges[5], oranges[9]))
helper = brewer.pal(12, "Paired")
col_fun_evidence = c("yes" = "#FFCC66", "no" = "white")
col_fun_coexpressed = c("yes" = "#FFCC66", "no" = "white")
col_fun_conserved = c("yes" = "#FFCC66", "no" = "white")
col_fun_coloc_gtex = c("yes" = "#FFCC66", "no" = "white")
col_fun_coloc_starnet = c("yes" = "#FFCC66", "no" = "white")
col_fun_bio = c("protein-coding" = "#FFC20A" , "non-coding RNA" = "#0C7BDC")

row_ha = rowAnnotation(foo4 = conserved, foo3 = coexpressed, foo5 = coloc_gtex, foo6 = coloc_starnet, foo = evidence, foo_empty = anno_empty(border = FALSE, width =unit(0.3, "cm")), foo2 = expression,# foo_bio = biotype,
	col = list(foo2 = col_fun_expression, foo = col_fun_evidence, foo3 = col_fun_coexpressed, foo4 = col_fun_conserved, foo5 = col_fun_coloc_gtex, foo6 = col_fun_coloc_starnet), # foo_bio = col_fun_bio),  ## set colors
	annotation_legend_param = list( ## define legend in detail
		foo = list(title = "supportive evidence", title_gp = gpar(fontsize = 14),labels_gp = gpar(fontsize = 12), legend_height = unit(2, "cm"), title_position = "topleft", ncol = 2, border = "grey"), 
		foo2 = list(title = "expression\nbreadth", title_gp = gpar(fontsize = 14),labels_gp = gpar(fontsize = 12), legend_width = unit(4, "cm"), title_position = "topleft", direction = "horizontal")), #,
	na_col = "white",
	simple_anno_size = unit(0.5, "cm"),
	show_legend = c(FALSE, FALSE, FALSE, FALSE,TRUE, FALSE, TRUE),#, TRUE),
	show_annotation_name = FALSE) ## remove annotation name 
	 
#border=FALSE needs to be set before specifying the axis_param
row_bar = rowAnnotation(foo_bio = biotype,  foo_empty = anno_empty(border = FALSE, width =unit(0.3, "cm")), bar2 = anno_barplot(mergedRems, gp = gpar(fill = "#CCCCCC", col ="grey"), bar_width = 0.85, border = FALSE, axis_param = list(gp=gpar(fontsize = 12))),width = unit(2, "cm"), show_annotation_name = FALSE, col = list(foo_bio = col_fun_bio), annotation_legend_param = list(foo_bio = list(title = "Biotype", title_gp = gpar(fontsize = 14),labels_gp = gpar(fontsize = 12), legend_height = unit(2, "cm"), title_position = "topleft",ncol = 2, boarder = "grey")))


ht_opt$ROW_ANNO_PADDING = unit(0.5, "cm")


ht = Heatmap(data_heatmap, name = "?",
	heatmap_legend_param = list(title = "number of CREs\nwith TF-SNVs", title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12), at = c(0,2 ,4,6,8), legend_width = unit(4, "cm"), title_position = "topleft", direction = "horizontal"),
	cluster_rows = FALSE, 
	cluster_columns = FALSE,
#	use_raster = TRUE, ## when switched on the white lines are gone but the matrix is blurry
	col = col_fun, 
	show_column_names = TRUE,
	show_column_dend  = FALSE,
	column_title = "cell type", 
	column_title_gp = gpar(fontsize = 20),
	column_names_gp = gpar(fontsize = 14),
	column_title_side = "bottom", #"top"
	show_row_names = TRUE,
	show_row_dend = FALSE,
	row_title = "Top 20 non-coding genes",
	row_title_gp = gpar(fontsize = 20),
	row_names_gp = gpar(fontsize = 14, fontface = "italic"),
	row_title_side = "left", #"right", #"left"
	left_annotation = row_ha,
	right_annotation = row_bar
)
draw(ht, heatmap_legend_side = "top")
decorate_annotation("bar2", {grid.text("CREs with TF-SNVs\nacross cell types", y = -0.2 ,x = 2, just = "bottom", gp=gpar(fontsize=14), rot = 0)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo", {grid.text("known CAD GWAS gene", y = -0.315 ,x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo2", {grid.text("expression breadth", y = -0.25 ,x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo3", {grid.text("coexpression", y = -0.17 , x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo5", {grid.text("colocalization GTEx", y = -0.255 , x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo6", {grid.text("colocalization STARNET", y = -0.305, x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo4", {grid.text("conserved in mouse", y = -0.255 , x = 0.85, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
decorate_annotation("foo_bio", {grid.text("Biotype", y = -0.1 , x = 0.75, just = "bottom", gp=gpar(fontsize=14), rot = 90)}) ## change label of the annotation via decorate_annotation
dev.off()

