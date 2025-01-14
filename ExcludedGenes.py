from collections import Counter
import HelperFunctions

"""Define a set of genes we want to exclude from all tests and analyses, namely pseudogenes, TEC genes, 
genes that don't have a unique 5'TSS and mirtrons (miRNAs that are entirely within gene bodies of other genes on 
the same strand)."""

annotation = "gencode.v38.annotation.gtf.gz"
out_file = "ExcludedGenes.txt"

excluded_genes = {}

# First get all the pseudogenes.
gene_biotypes = HelperFunctions.gene_biotypes(annotation)
for gene, val in gene_biotypes.items():
    if val['general'] == 'pseudogene':
        excluded_genes[gene] = {'pseudogene'}
    if val['general'] == 'TEC':
        excluded_genes[gene] = {'TEC'}

# Secondly, find all genes that don't have a unique 5'TSS.
gene_tss = HelperFunctions.gene_window_bed(annotation, extend=200, gene_set=set(), tss_type='5', dict_only=True)
tss_counter = Counter([val['chr'] + '-' + str(next(iter(val['tss']))) for val in gene_tss.values()])
non_unique_genes = set([g for g, val in gene_tss.items() if tss_counter[val['chr'] + '-' + str(next(iter(val['tss'])))] > 1])
for gene in non_unique_genes:
    if gene not in excluded_genes:
        excluded_genes[gene] = {'non-unique 5TSS'}
    else:
        excluded_genes[gene].add('non-unique 5TSS')

# Check if we can find mirtrons, meaning they are located in other genes. We don't restrict to introns though.
mirnas = set([g for g, vals in gene_biotypes.items() if vals['gtf'] == 'miRNA'])
mirna_gbs = HelperFunctions.gene_body_bed(annotation, gene_set=mirnas, dict_only=False)
all_nonmirna_gbs = HelperFunctions.gene_body_bed(annotation, gene_set=set([g for g in gene_biotypes.keys() if g not in mirnas]))
mirna_in_gbs = set([x.fields[3] for x in mirna_gbs.intersect(all_nonmirna_gbs, f=1, s=True, u=True)])
for gene in mirna_in_gbs:
    if gene not in excluded_genes:
        excluded_genes[gene] = {'miRNA within host gene'}
    else:
        excluded_genes[gene].add('miRNA within host gene')

with open(out_file, 'w') as gene_out:
    gene_out.write('\n'.join([g + '\t' + ','.join(val) for g, val in excluded_genes.items()]))


