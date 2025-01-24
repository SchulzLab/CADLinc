import statsmodels.stats.multitest
import pandas as pd
import pickle

"""Read the GATES results, filter genes, do FDR-adjustment and set conditions for genes to be considered relevant."""

open_genes_pkl = "OpenGenesCell.pkl"
open_genes = pickle.load(open(open_genes_pkl, 'rb'))
all_open = set().union(*open_genes.values())
excluded_genes_file = "ExcludedGenes.txt"
excluded_genes = set([x.strip().split('\t')[0] for x in open(excluded_genes_file).readlines()])

# Checking FDR for GATES.
gates_file = "result_gene_based_test_genebody_Filtered.txt"
gates_df = pd.read_table(gates_file, sep=',', header=0)
gates_df = gates_df[(gates_df['gene_list...1.'].isin(all_open)) & (~gates_df['gene_list...1.'].isin(excluded_genes))]
own_fdr = statsmodels.stats.multitest.multipletests(gates_df['gates_df'].values, alpha=0.05, method='fdr_bh')[1]
gates_df['FDR'] = own_fdr
these_genes = set(gates_df[(gates_df['FDR'] <= 0.025) & (gates_df['no_snps_per_gene'] >= 10)]['gene_list...1.'])
gates_df['Filtered'] = [g in these_genes for g in gates_df['gene_list...1.']]
gates_df.to_csv(gates_file.replace(".txt", '_Filtered.txt'), sep='\t', header=True, index=False)



