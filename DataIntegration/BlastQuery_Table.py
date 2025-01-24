# https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html

"""Format the blast output into one table."""

# Show plots as part of the notebook and make tools available
# %matplotlib inline
import matplotlib
matplotlib.use('Agg')
matplotlib.style.use('ggplot')

# Standard library packages
import os

# Import Numpy, Pandas and Seaborn
import pandas as pd

fwd_out = os.path.join('mm_alignments_lncRNA.txt')
rev_out = os.path.join('hg_alignments_lncRNA.txt')

# Load the BLAST results into Pandas dataframes
fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
rev_results = pd.read_csv(rev_out, sep="\t", header=None)

# Add headers to forward and reverse results dataframes
headers = ["query", "subject", "identity", "coverage",
           "qlength", "slength", "alength",
           "bitscore", "E-value"]
fwd_results.columns = headers
rev_results.columns = headers

# Create a new column in both dataframes: normalised bitscore
fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

# Create query and subject coverage columns in both dataframes
fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
rev_results['qcov'] = rev_results.alength/rev_results.qlength
fwd_results['scov'] = fwd_results.alength/fwd_results.slength
rev_results['scov'] = rev_results.alength/rev_results.slength

# Clip maximum coverage values at 1.0
fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
rev_results['scov'] = rev_results['scov'].clip(upper=1)

# Inspect the forward results data
print("Forward alignment Results ...........\n",fwd_results.head())

# Inspect the reverse results data
print("Reverse alignment Results ...........\n",rev_results.head())

#=========================== for gene level: forward alignments ===================================#

mouse_anno = os.path.join('annotations', 'mm_lncRNA_transcript_gene.txt')
mouse_anno = pd.read_csv(mouse_anno, sep="\t", header=0)
mouse_anno.drop(columns=['transcript_biotype', 'external_synonym'], inplace=True)

# print("Moues Annotations-------------------\n",mouse_anno.head())

fwd_results_mousegene = fwd_results.merge(mouse_anno.rename({'ensembl_transcript_id_version': 'mouse_query', 
                                          'external_gene_name':'Mouse_query_lncgene', # 'external_synonym':'Mouse_synonym',
                                          'ensembl_gene_id_version': 'query_gene'}, axis=1), 
                left_on='query', 
                right_on='mouse_query',
                how='left')
fwd_results_mousegene = fwd_results_mousegene.drop_duplicates()

human_anno = os.path.join('annotations', 'hg_lncRNA_transcript_gene.txt')
human_anno = pd.read_csv(human_anno, sep="\t", header=0)
human_anno.drop(columns=['transcript_biotype', 'ensembl_transcript_id', 'external_synonym'], inplace=True)


fwd_results_msgn_hmgn = fwd_results_mousegene.merge(human_anno.rename({
                                                    'ensembl_transcript_id_version': 'human_subject', 
                                                    'external_gene_name':'human_lncgene',  # 'external_synonym':'Human_synonym',
                                                    'ensembl_gene_id_version': 'subject_gene'}, axis=1),
                                        left_on='subject', 
                                        right_on='human_subject',
                                        how='left')
fwd_results_msgn_hmgn = fwd_results_msgn_hmgn.drop_duplicates()


# ============================= for gene level: reverse alignments ===============================

rev_results_mousegene = rev_results.merge(mouse_anno.rename({
                                            'ensembl_transcript_id_version': 'mouse_subject', 
                                            'external_gene_name':'mouse_query_lncgene', # 'external_synonym':'Mouse_synonym',
                                            'ensembl_gene_id_version': 'subject_gene'}, axis=1), 
                left_on='subject', 
                right_on='mouse_subject',
                how='left')

rev_results_mousegene = rev_results_mousegene.drop_duplicates()

rev_results_msgn_hmgn = rev_results_mousegene.merge(human_anno.rename({
                                            'ensembl_transcript_id_version': 'human_query', 
                                            'external_gene_name':'human_query_lncgene', # 'external_synonym':'Mouse_synonym',
                                            'ensembl_gene_id_version': 'query_gene'}, axis=1), 
                left_on='query', 
                right_on='human_query',
                how='left')
rev_results_msgn_hmgn = rev_results_msgn_hmgn.drop_duplicates()
print("Merged reverse alignments with mouse and human annotations .......\n",rev_results_msgn_hmgn.head())

#===============================================================================================#
# Merge forward and reverse results
rbbh_gene = pd.merge(fwd_results_msgn_hmgn, rev_results_msgn_hmgn[['query_gene', 'subject_gene']],
                left_on='subject_gene', right_on='query_gene',
                how='outer')
rbbh_gene.to_csv('rbbh_gene_intermediate.csv', index=False, sep="\t")
# Discard rows that are not RBH
rbbh_gene = rbbh_gene.loc[rbbh_gene.query_gene_x == rbbh_gene.subject_gene_y]
rbbh_gene = rbbh_gene.drop_duplicates()
rbbh_gene.to_csv('rbbh_gene_intermediateII.csv', index=False, sep="\t")


# Group duplicate RBH rows, taking the maximum value in each column
rbbh_gene = rbbh_gene.groupby(['query_gene_x', 'mouse_query', 'human_subject', 'subject_gene_x', ]).max()

# Inspect the results
print("========================== RBBH gene results ========================\n",rbbh_gene.head())

rbbh_gene.to_csv('rbbh_gene_alignments.csv', index=False, sep="\t")
