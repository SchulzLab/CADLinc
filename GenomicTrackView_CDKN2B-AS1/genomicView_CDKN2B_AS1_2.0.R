library(Gviz)
library(GenomicFeatures)

library(org.Hs.eg.db)

##hg38

#input files
args <- commandArgs(trailingOnly = TRUE)
inputFile <- args[1] # bed file with cell type specific interactions IQCH-AS1
snvFile <- args[2]
annotationFile <- args[3] 
output <- args[4]

pdf(output) ## open output file

# genomic regions our samples
data_ <- data.frame(read.delim(inputFile, header = TRUE, sep = "\t")) # read in bed file (needs the header chr start end strand id usw)
print(head(data_))
data_SNVs <- data.frame(read.delim(snvFile, header = TRUE, sep = "\t")) # read in bed file (needs the header chr start end strand id usw)
print(head(data_SNVs))

chr <- "chr9"

atrack = AnnotationTrack(range = data_, name = "regulatory elements (CREs)",chromosome = chr,  genome = "hg38", col.line = "white", col = "white", ldw = 0.5, shape = "box") ## annotations for IQCH-AS1
feature(atrack) <- ifelse(strand(atrack)=="+", "plus", "minus") 

## SNVs 
track_SNV = AnnotationTrack(range = data_SNVs, name = "regulatory SNVs",chromosome = chr,  genome = "hg38", col.line = "white", col = "white", ldw = 0.5, shape = "box")

## only lead snvs in region 2195000 to 22150000
feature(track_SNV) <- c( "rs1333039", "rs1537373",  "rs3731238",  "rs4977753", "rs61271866", "rs76959412")

gtrack <- GenomeAxisTrack() ##genomic trace

geneModels = loadDb("gencode.v38.annotation.TxDb")
print(columns(geneModels))

grtrack_IQCH_AS1 <- GeneRegionTrack(geneModels,genome = "hg38", chromosome = "chr9", name = "Gene annotation", collapseTranscripts = "meta", shape = "smallArrow", transcriptAnnotation = "symbol", fill = "#0C7BDC", ldw = 0.5, col = "#0C7BDC", col.line = "#0C7BDC")
print("here")
print(grtrack_IQCH_AS1)

## get the gene names 
z <- ranges(grtrack_IQCH_AS1)
print(head(z$gene))
z$symbol <- mapIds(x = org.Hs.eg.db, keys =gsub("\\..*","",z$gene), keytype = "ENSEMBL",column = "SYMBOL") ## has to remove the version number after the dot which is done with gsub(...)
z$symbol <- ifelse(is.na(z$symbol), z$gene, z$symbol)
print(head(z))
ranges(grtrack_IQCH_AS1) <- z

## snvs right side 
plotTracks(list(track_SNV, atrack,  grtrack_IQCH_AS1, gtrack),  from = 21970000 , to = 22110000 , plus="#DC3220", minus="azure4", rs2165408 = "black", rs12352425 = "black", rs1333039 = "black", rs1537373 = "black", rs3217986 = "black", rs3731238 = "black", rs3928893 = "black", rs4977753 = "black", rs61271866 = "black", rs76959412 = "black", rs77777666 = "black"  ,groupAnnotation = "group", sizes = c(0.2,6.0,1, 0.5))

dev.off()
