library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

#log2FC_threhsold=0.3

tissues = c("AOR", "LIV", "SF", "SKLM", "VAF")

sigTFs = data.frame()
for (t in tissues) {

	print(t)
	dat <- read.csv(paste("TFs_DEGs_", t, "_extended.csv", sep = ""), header = TRUE)


	# Create a new column "delabel2" to de, that will contain the name of enriched TFs (NA in case they are not)
	dat$delabel2 <- NA
	dat$delabel2[dat$enrichedTF != "FALSE"] <- dat$Gene.name[dat$enrichedTF != "FALSE"]

	# add a column of NAs
	dat$diffexpressed2 <- "NO"
	dat$diffexpressed2[dat$padj <= 0.05] <- "YES"

	### compute hypergeometric test 
	## using the R function phyper

	## collect the subsets 
	# q: number of white balls drawan without replacement from an urn which contains both white and black balls -> number of TFs with high odds-ratio that are DEG
	TFsHighOddsRatioDEG = subset(dat, dat[, "diffexpressed2"] == "YES" & dat[,"enrichedTF"] == "TRUE")
	TFsHighOddsRatioDEG$tissue <- t
	sigTFs = rbind(TFsHighOddsRatioDEG, sigTFs)
	print(nrow(sigTFs))
	
	print(paste("TFs with high Odds ratio DEGs: ",nrow(TFsHighOddsRatioDEG), sep = ""))

	# m: number of white balls in the urn -> all TFs with high-odds ratio -> 34
	TFsHighOddsRatio = subset(dat, dat[,"enrichedTF"] == "TRUE")
	print(paste("TFs high odds ratio: ", nrow(TFsHighOddsRatio), sep = ""))

	# n: number of black balls in the urn -> number of TFs NO high-odds ratio -> 
	TFsNoHighOddsRatio = subset(dat, dat[,"enrichedTF"] == "FALSE")

	print(paste("TFs no high odds ratio: ", nrow(TFsNoHighOddsRatio), sep = ""))

	# k: number of balls drawn from the urn -> Number of TFs that are DEG
	TFsDEG =  subset(dat, dat[, "diffexpressed2"] == "YES")
	print(paste("TFs DEG: " , nrow(TFsDEG), sep = ""))

	##perfrom the test: 

	hypergeo_test=phyper(q =nrow(TFsHighOddsRatioDEG) , m = nrow(TFsHighOddsRatio) , n= nrow(TFsNoHighOddsRatio) , k = nrow(TFsDEG) , lower.tail=FALSE, log.p = FALSE)
	print(hypergeo_test)

	## add color scale 
	mycolors <- c("grey", "orange")

	ggplot(data=dat, aes(x=log2FoldChange, y=-log10(padj), col = enrichedTF, label = delabel2)) + 
		geom_point() + 
		geom_text_repel() + 
    		geom_hline(yintercept=-log10(0.05), col="black", linetype="dotted") + ## add horizontal line
		scale_colour_manual(values = mycolors) + 
		ggtitle(paste("Tissue: ",t, ", p-value hypergeometric test: ", round(hypergeo_test,digits = 4), sep = "")) + 
		labs(col = "CAD-TF") + 
		xlab("log2FoldChange") + 
		ylab("-log10(FDR)") + 
		theme_minimal() + 
		theme(text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) 

	ggsave(paste("vulcanoPlot_TFs_",t,".pdf", sep = ""))
}
print("sigTFs")
print(head(sigTFs))

## collecting all sigTFs over all tissues and plot them as piechart 
## 34 is the number of TFs with high odds ratio

tfs = unique(sigTFs[,"delabel2"])
print("piechart")
print(length(tfs))

pie_data = data.frame(group = c("non-differentially expressed CAD-TFs","differentially expressed CAD-TFs"), value = c(34 - length(tfs),length(tfs)))
ggplot(pie_data, aes(x="", y=value, fill=group)) +
  	geom_bar(stat="identity", width=1) +
	scale_fill_manual(name = "CAD-TF expression", labels = c("non-differentially expressed CAD-TFs" = "not differential (2/34)", "differentially expressed CAD-TFs" =  "differential (32/34)"), values=c("#77AADD","#999999")) + 
  	coord_polar("y", start=0) + 
	theme_void() + 
	theme( legend.title=element_text(size=20), legend.text=element_text(size=18))

ggsave("piechart_allSTARNET_tissue.pdf")


counts_ = data.frame(table(sigTFs$tissue))
counts_$percentage = counts_$Freq/34

sigTFs$diffexpressed3 <- "NO"
sigTFs$diffexpressed3[sigTFs$padj < 0.01 & sigTFs$log2FoldChange > 0.3] <- "UP"
sigTFs$diffexpressed3[sigTFs$padj < 0.01 & sigTFs$log2FoldChange < -0.3] <- "DOWN"

counts_overall = data.frame()
for (t in tissues) {
	tissue_specific =  subset(sigTFs, sigTFs[, "tissue"] == t)
	counts = data.frame(table(tissue_specific$diffexpressed3))
	helper = data.frame(tissue = c(t,t,t), diff = c("1Up", "2down", "3none"), count = c(counts[2,"Freq"], counts[1,"Freq"], 39 - counts[2,"Freq"] - counts[1,"Freq"]))
	counts_overall = rbind(helper, counts_overall)
}
print("counts overall")
print(counts_overall)
write_delim(counts_overall, "sigTFs_perTissue_26_07.txt")
## barplot per celltype

ggplot(data=counts_overall, aes(x=reorder(tissue, count), y=count, fill = diff)) +
  	geom_bar(stat="identity", alpha = 0.8)+
	xlab("STARNET tissue") + 
	ylab("Number of CAD-TFs") + 
	ylim(0,40) + 
	scale_x_discrete(labels=c("AOR" = 'atherosclerotic\naortic root', "LIV" = 'liver', "SF" = 'subcutaneous fat', "SKLM" = 'skeletal muscle', "VAF" = 'visceral fat'), guide = guide_axis(angle = 90)) + 
	scale_fill_manual("gene expression", labels = c("1Up" = "up", "2down" = "down", "3none"= "not differential"), values = c("#228833", "#CCBB44", "#999999")) +
	theme_classic() + 
	theme(text = element_text(size = 18), axis.text.x = element_text(hjust = 1, size = 16), legend.position="top") 

ggsave("barplot_STARNET_tissue_up_down.pdf")

ggplot(data=counts_, aes(x=reorder(Var1, Freq), y=Freq)) +
  	geom_bar(stat="identity", fill = "grey", alpha = 0.8)+
	xlab("STARNET tissue") + 
	ylab("Number of differentially expressed CAD-TFs") + 
	ylim(0,34) + 
	scale_x_discrete(labels=c("AOR" = 'atherosclerotic aortic root', "LIV" = 'liver', "SF" = 'subcutaneous fat', "SKLM" = 'skeletal muscle', "VAF" = 'visceral fat'), guide = guide_axis(angle = 45)) + 
	theme_minimal() + 
	theme(text = element_text(size = 14), axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) 

ggsave("barplot_STARNET_tissue.pdf")


sigTFs = sigTFs[c("Ensembl_ID", "baseMean", "log2FoldChange",  "lfcSE","stat", "pvalue", "padj", "Gene.name", "X.log10.pvalue_adj.", "enrichedTF", "diffexpressed2", "diffexpressed3", "tissue")]
print(head(sigTFs))

write_delim(sigTFs, "sigTFs_perTissue.txt")














