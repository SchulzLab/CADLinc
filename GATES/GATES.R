#install.packages("COMBAT")

library(COMBAT)
library(dplyr)

# Function of GATES
gates_f = function(x, cor_G){

  pval_sort  <- sort(x)
  pval_order <- order(x)
  n_snps     <- length(x)

  cor_P <- cor_G[pval_order, pval_order]
  cor_P <- 0.2982*cor_P^6 - 0.0127*cor_P^5 + 0.0588*cor_P^4 + 0.0099*cor_P^3 + 0.6281*cor_P^2 - 0.0009*cor_P
  
  p_gates <- ext_simes_f(pval_sort, cor_P)

  p_gates
  
}

# Function of extended Simes
ext_simes_f = function(x, cor_r){
  eff.snpcount.fun <- function(ldmat) {
    ldmat <- as.matrix(ldmat)
    snpcount.local <- dim(ldmat)[1]
    if (snpcount.local <= 1) return(1)
    ev <- eigen(ldmat, only.values = TRUE)$values
    if (sum(ev < 0) != 0) {
      ev <- ev[ev > 0]
      ev <- ev/sum(ev) * snpcount.local
    }
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  eff.snpcount.global <- eff.snpcount.fun(cor_r)
  
  n_values <- length(x)
  candid <- sapply(1:n_values, function(i){
    (eff.snpcount.global * x[i])/eff.snpcount.fun(cor_r[1:i,1:i])
    
  })
  
  p_ext_simes <- min(candid)
  p_ext_simes
}

### INPUT
args <- commandArgs(trailingOnly = TRUE)
inputSNPs <- args[1] 
consideredGenes <- args[2] 
LD_structure <- args[3]
outputFile <- args[4]

#### Read P-values and SNPs
data_snp_p <- read.csv2(inputSNPs, header=TRUE,sep="\t")

snp.info   <- as.data.frame(data_snp_p)
snp.info <- snp.info[c('SNP_ID','P.value','GENE')]

#### Read gene list
gene_list <- read.csv2(consideredGenes, header=TRUE,sep="\n") ###Pablos genes of interest
gene.info   <- as.data.frame(gene_list[,1])

### For all genes
snp.pvals_ <- as.data.frame(snp.info)

#### Initialize dataframes 
gates_df <- matrix(ncol=1, nrow=nrow(gene.info))
no_snps_per_gene <- matrix(ncol=1, nrow=nrow(gene.info))

for (i in 1:nrow(gene.info))
{
  gene_to_test <- gene.info[[1]][i]
  snp.pvals <- snp.pvals_[grepl(gene_to_test, snp.pvals_$GENE),]
  ### check if gene has important SNPs
  if (dim(snp.pvals)[1] != 0)
  {
    snplist <- list(snp.pvals$SNP_ID)
    
    snp.pvals <- as.matrix(snp.pvals)
    snp.pvals<- subset(snp.pvals, select = - c(GENE))
    snp.pvals <- subset(snp.pvals, select = - c(SNP_ID))
    snp.pvals <- as.numeric(snp.pvals)

    if (length(snplist[[1]]) == 1)
    {
      pval_gates <- snp.pvals
      
    }
    else
    {
    ###ignore genes with less than 3 SNPs
    gene_LD <- paste(gene_to_test,"txt","ld",sep=".")
    ### Read LD R SQUARED
    path_to_LD <- paste(LD_structure,gene_LD,sep="")

    if (file.exists(path_to_LD)){
	print("path exists")
	print(path_to_LD)
    	cor_G <- read.csv2(path_to_LD, header = FALSE,sep=" ")
    	cor_G   <- as.matrix(cor_G)
    	cor_G   <- as.data.frame(cor_G)
    	cor_G <- cor_G[1:(length(cor_G)-1)]
    	cor_G <- data.frame(sapply(cor_G, function(x) as.numeric(as.character(x))))

    	if (length(snplist[[1]]) > dim(cor_G)[1])
    	{
    	  pval_gates <- 888888
      
    	}
    	else
    	{
        	 #### GATES TEST
    	pval_gates <- gates_f(x=snp.pvals, cor_G=cor_G)
    	}
    }
	else
	{
		print("no LD structure for gene")
		print(path_to_LD)
    }
}
    ##gate df
    number_of_snps <- length(snplist[[1]])
    gates_df[i,] <- pval_gates
      ### number of snps
    no_snps_per_gene[i,] <- number_of_snps

  }  
}

### Concatenate info
result_data <- data.frame(gene.info,gates_df,no_snps_per_gene)

result_data <- result_data[order(result_data$gates_df),]
## multiple testing corrected p-values
p.adjuste_values_G <- p.adjust(result_data$gates_df, method = p.adjust.methods, n = length(result_data$gates_df))
result_data$GATES_adj <- p.adjuste_values_G

write.csv(x=result_data, file=outputFile)



