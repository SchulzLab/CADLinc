load("sneep_input.RData") ## Loads the variables: "genes", "x", "x_annotation"
print("read first file")
load("coexpression_PearsonCor.RData") # Loads the variable: "C_res"
print("read second file")

### The x contains the gene expression of 9664 samples from GTEx
### The C_res contains the pairwise correlation between the genes, but it does not have the row/colnames. So, we need the variable "x" from about to form the connection
top <- 1000
genes <- readLines("250924_JointKnownGenes_EnsemblIDs.txt")
res <- list();
for(gn in genes){
	 hit <- which(x$gene_id == gn);
	if(length(hit)){
		res[[gn]] <- x$gene_id[order(C_res[hit, ], decreasing=T)] ## Take the genes that are most correlated with the gene of interest stored in genes
	}
}
res_mat <- matrix(NA, ncol= top + 1, nrow= length(res))
for(i in seq(length(res))){
	res_mat[i, ] <- c(names(res)[i], res[[i]][seq(2, top+1)])
}
write.table(res_mat, paste0("co_expressed_genes_JointKnownGenes_2509_", top, ".txt"), sep= "\t", quote= F, row.names= F, col.names= F) ## Write the results in the file

