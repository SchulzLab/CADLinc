import pickle
import pandas as pd
from itertools import chain
import HelperFunctions


"""Run GO enrichment on the target genes of genes."""

plot_tag = ''
plot_out = "Plots/" + plot_tag + '_'

excluded_genes_file = "ExcludedGenes.txt"
excluded_genes = set([x.strip().split('\t')[0] for x in open(excluded_genes_file).readlines()])

open_genes_pkl = "OpenGenesCell.pkl"
open_genes = pickle.load(open(open_genes_pkl, 'rb'))
all_open_genes = set(chain(*open_genes.values())) - excluded_genes

tf_file = "TF_table_2905_extended.txt"
tf_table = pd.read_table(tf_file, sep='\t', header=0)

open_tfs = tf_table[tf_table['openPromoter']]
selected_tfs = open_tfs['TF'].values
all_tfs_targets = {}

for tf in selected_tfs:
    print(tf)
    table_hit = open_tfs[open_tfs['TF'] == tf]
    tf_targets = set(table_hit['targetGenes'].values[0].strip(',').split(',')) - excluded_genes
    target_dict = {tf: tf_targets}
    all_tfs_targets[tf] = tf_targets

# Updated version with corrected target gene sets.
tf_groups = [{'tfs': ['PATZ1', 'CTCF(MA0139.1)', 'CTCFL', 'NFIX(MA1528.1)', 'TFAP2E'], 'sources': ['GO:BP'], 'out_tag': 'Immune response and inflammation'},
             {'tfs': ['TCFL5'], 'sources': ['GO:MF'], 'out_tag': 'Metabolite transportation'},
             {'tfs': ['TCFL5'], 'sources': ['GO:BP'], 'out_tag': 'Innervation and synaptic transmission'},
             {'tfs': ['KLF16'], 'sources': ['WP'], 'out_tag': 'TGF beta signaling pathway'},
             {'tfs': ['THRB(MA1575.1)', 'PATZ1'], 'sources': ['REAC'], 'out_tag': 'Lipid metabolism'},
             {'tfs': ['THRB(MA1575.1)', 'TFEB'], 'sources': ['WP'], 'out_tag': 'Lipid metabolism'},
             {'tfs': ['THRB(MA1575.1)'], 'sources': ['GO:BP'], 'out_tag': 'Lipid metabolism'},
             {'tfs': ['ZNF610'], 'sources': ['REAC'], 'out_tag': 'Vascular remodeling'}]
for group in tf_groups:
    # Remove motif ID.
    group_tf_targets = {t.split('(')[0]: all_tfs_targets[t] for t in group['tfs']}
    _ = HelperFunctions.go_enrichment(group_tf_targets, title_tag=group['out_tag'], numerate=True,
                                   out_tag=plot_out+group['out_tag'], max_terms='all', wanted_sources=group['sources'],
                                   background={t: all_open_genes for t in group_tf_targets})
