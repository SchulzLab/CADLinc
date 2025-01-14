library(biomaRt)
# choosing the dataset
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "asia")

filters = listFilters(ensembl)
searchFilters(mart = ensembl, pattern = "transcript_biotype")

# ========================= (incase) obj: get the lncRNA ids : for Human =============================
lnc_chr = getBM(attributes = c('entrezgene_id','transcript_biotype','ensembl_gene_id', 'hgnc_symbol'),
      filters = c('chromosome_name', 'transcript_biotype'),
      values = list("lncRNA"),
      mart = ensembl)

# get the sequences of the entrez ids: for Human
lncRNAseqchr = getSequence(id=lnc_chr$hgnc_symbol, #entrezgene_id,
                              type="hgnc_symbol",
                              seqType="gene_exon_intron",
                              mart=ensembl)

write.table(lncRNAseqchr, "ensembl/annotation/Homo_sapiens.GRCh38.single_line.lncrna.txt", quote=F)
exportFASTA(lncRNAseqchr, "ensembl/annotation/Homo_sapiens.GRCh38.single_line.lncrna.fa")

# ================== obj: get the lncRNA ids : for Mouse ====================
listMarts()
ensemblmm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

mm_lnc_chr = getBM(attributes = c('entrezgene_id','transcript_biotype','ensembl_gene_id', 'external_gene_name'),
                  filters = c('transcript_biotype'),
                  values = list("lncRNA"),
                  mart = ensemblmm)

# get the sequences of the entrez ids: for Mouse
mm_lncRNAseqchr = getSequence(id=mm_lnc_chr$external_gene_name,
                              type="external_gene_name",
                             seqType="gene_exon_intron",
                             mart=ensemblmm)

write.table(mm_lncRNAseqchr, "ensembl/annotation/Mus_musculus.GRCm39.single_line.lncrna.txt", quote=F)
exportFASTA(mm_lncRNAseqchr, "ensembl/annotation/Mus_musculus.GRCm39.single_line.lncrna.fa")


# ============================== algorithm RBH =============================#
# lncRNA sequences downloaded from GENCODE

system("makeblastdb -in Mus_musculus.GRCm39.single_line.lncrna.fa -out mouse_blastidx_lncRNA -parse_seqids -dbtype nucl")

system("makeblastdb -in Homo_sapiens.GRCh38.single_line.lncrna.fa -out human_blastidx_lncRNA -parse_seqids -dbtype nucl")


system("blastn -query Mus_musculus.GRCm39.single_line.lncrna.fa -db human_blastidx_lncRNA -num_threads 20 -outfmt \"6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\" -max_target_seqs 1 -out mm_alignments_lncRNA.txt")


system("blastn -query Homo_sapiens.GRCh38.single_line.lncrna.fa -db mouse_blastidx_lncRNA -num_threads 20 -outfmt \"6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\" -max_target_seqs 1 -out hg_alignments_lncRNA.txt")



