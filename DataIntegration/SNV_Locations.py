from pybedtools import BedTool
import PyComplexHeatmap as pch
from PyComplexHeatmap import *
import matplotlib.pylab as plt
import HelperFunctions
import DataIntegration.PlottingFunctions as PlottingFunctions

"""Plot pie charts of the location of SNPs in relation to genes, and also the number of SNPs in
enhancers of each cell.
CARE: Don't use bedtool's merge on the SNPs, it will also merge neighbouring SNPs."""

plot_out = "/Plots/"
annotation = "gencode.v38.annotation.gtf.gz"
enhancer_folder = "enhancer_interactions_avghg38_unique/"
gene_table_file = "2509_JointKnownGenes_GeneTable.txt"

all_rems_file = "AllREMsMerged_200624.bed"
if not os.path.isfile(all_rems_file):
    rems = []
    for file in [x for x in os.listdir(enhancer_folder) if x.endswith(".bed")]:
        for entry in open(enhancer_folder + '/' + file).readlines():
            rems.append('\t'.join(entry.strip().split('\t')[:3] + [file.split('_')[0]]))
    merged_rems = BedTool('\n'.join(rems), from_string=True).sort().merge(c=[4], o=['distinct'])
    open(all_rems_file, 'w').write(str(merged_rems))

# Get a bed-object for each SNP set.
snp_beds = {}
# First all significant ones from the GWAS, threshold as used for SNEEP, including the proxy SNPs.
ld_file = "GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_allLDSNVs.txt"
ld_df = pd.read_table(ld_file, header=0, sep='\t')
snp_beds['1%FDR SNVs & proxies'] = BedTool('\n'.join(set(['\t'.join([str(y) for y in x]) for x in ld_df[['chr', 'start', 'end']].values])), from_string=True)

# Get the regulatory SNPs from SNEEP.
sneep_file = "sneep_cadlincGWAS_2.52x0.00001_22_05.txt"
sneep_df = pd.read_table(sneep_file, header=0, sep='\t')
snp_beds['TF-SNVs'] = BedTool('\n'.join(set(['\t'.join([x.split(':')[0], x.split(':')[1].split('-')[0], x.split('-')[1]])
                                        for x in sneep_df['SNP_position'].values])), from_string=True)
snp_beds['TF-SNVs in CREs'] = snp_beds['TF-SNVs'].intersect(all_rems_file, u=True)
cadsnvs_in_cres = snp_beds['1%FDR SNVs & proxies'].intersect(BedTool(all_rems_file), u=True)
snp_sets = {k: set([str(x).strip() for x in vals]) for k, vals in snp_beds.items()}

# A Venn of the CAD-SNVs, the TF-SNVs and which are in CREs.
venn_snps = {'CAD-SNVs': snp_sets['1%FDR SNVs & proxies'], 'TF-SNVs': snp_sets['TF-SNVs'], 'in CREs': set([str(x).strip() for x in cadsnvs_in_cres])}
PlottingFunctions.venn_from_list(list(venn_snps.values()), ['' for k in venn_snps], plot_path=plot_out+"SNPsVenn",
                            title='', scaled=True, number_size=12, xsize=10, ysize=5)

snp_locs = HelperFunctions.gene_location_bpwise(gtf_file=annotation, bed_dict=snp_beds, plot_path=plot_out, tss_type='5')

# ---------------------------------------------------------------------------------------------------------------------
# Barplot intersection SNPs and enhancers
# ---------------------------------------------------------------------------------------------------------------------
# Plot a barplot with the number of SNPs in enhancers of each cell type.
meta_file = "Enhancer_metadata.txt"
meta_df = pd.read_table(meta_file, sep='\t', header=0)
cell_types = meta_df['ID'].values
cell_name_map = {val['ID']: val['Label'] for val in meta_df.to_dict(orient='records')}

enhancer_snps = []
for cell in cell_types:
    print(cell)
    match_file = [enhancer_folder + '/' + x for x in os.listdir(enhancer_folder) if x.startswith(cell.split('_')[0]+'_') and 'hg38' in x][0]
    these_enhancer = BedTool(match_file).sort().merge()
    for snp_t, snp_b in snp_beds.items():
        enh_w_snp = len(these_enhancer.intersect(snp_b, u=True))
        enhancer_snps.append([cell, snp_t, len(snp_b.intersect(these_enhancer, u=True)), enh_w_snp, enh_w_snp / len(these_enhancer)])

reorder_dict = {c: {} for c in cell_types}
for entry in enhancer_snps:
    reorder_dict[entry[0]][entry[1]] = entry[4]
    reorder_dict[entry[0]][entry[1] + ' in CREs'] = entry[2]  # NOTE here the TF-SNVs in CREs is getting overwritten.

cell_enhancer_snp_df = pd.DataFrame.from_dict(reorder_dict, orient='index')
cell_enhancer_snp_df['Cell type'] = [cell_name_map[i] for i in cell_enhancer_snp_df.index]
cell_enhancer_snp_df = cell_enhancer_snp_df.join(meta_df.set_index('ID')[['Class']])

# Stacked barplots with the absolute number of 1%FDR SNVs in REMs and how many of those are regSNVs. Need to reformat
# the df again as the regSNVs are a subset of the 1%FDR.
cell_enhancer_snp_forstack = cell_enhancer_snp_df[['1%FDR SNVs & proxies in CREs', 'TF-SNVs in CREs in CREs']]
cell_enhancer_snp_forstack.columns = ['1%FDR SNVs & proxies', 'TF-SNVs']
cell_enhancer_snp_forstack['non TF-SNVs'] = cell_enhancer_snp_forstack['1%FDR SNVs & proxies'] - cell_enhancer_snp_forstack['TF-SNVs']
cell_enhancer_snp_forstack['Cell type'] = [cell_name_map[i] for i in cell_enhancer_snp_forstack.index]
cell_enhancer_snp_df['rel TF-SNVs'] = cell_enhancer_snp_df['TF-SNVs']
joint_cell_df = cell_enhancer_snp_df[['1%FDR SNVs & proxies', 'rel TF-SNVs', 'Class']].join(cell_enhancer_snp_forstack[['TF-SNVs', 'non TF-SNVs']])

# Sort lineages by highest relative fraction.
lineage_order = list(joint_cell_df.groupby('Class')[['1%FDR SNVs & proxies']].mean().sort_values(by='1%FDR SNVs & proxies', ascending=True).index)
joint_cell_df['lin_index'] = [lineage_order.index(l) for l in joint_cell_df['Class']]
joint_cell_df = joint_cell_df.sort_values(by=["lin_index", '1%FDR SNVs & proxies'], ascending=False)
joint_cell_df.index = [cell_name_map[c] for c in joint_cell_df.index]

all_groups = list(dict.fromkeys(joint_cell_df['Class']).keys())

# This part is admittedly very ugly, as I never managed to convince PyComplexHeatmap to really plot what I wanted, and
# instead had to plot parts of it which I later puzzled together.
colour_dict = {c: PlottingFunctions.glasbey_palettes['glasbey_cool'][i] for i, c in enumerate(all_groups)}
left_ha = pch.HeatmapAnnotation(label=pch.anno_label(joint_cell_df['Class'], merge=True, legend=False, extend=False, colors='black', relpos=(1, 0.5)),
                                Group=pch.anno_simple(joint_cell_df['Class'], legend=False, colors=colour_dict),
                                verbose=1, axis=0, plot_legend=False, label_kws=dict(visible=False))
right_annotation = HeatmapAnnotation(axis=0,orientation='right',
                                Row=anno_barplot(joint_cell_df[['TF-SNVs', 'non TF-SNVs']], colors=['#C52A1C', '#c7c7c7'],legend=True, height=60, linewidth=0.1),
                                verbose=0, label_side='top', label_kws={'horizontalalignment': 'left','rotation':45,'visible':False})
plt.figure(figsize=(5, 10))
cm = pch.ClusterMapPlotter(data=joint_cell_df[['1%FDR SNVs & proxies', 'rel TF-SNVs']],# left_annotation=left_ha,
                           right_annotation=right_annotation, row_cluster=False, col_cluster=False,
                       label='Fraction of SNV type', row_dendrogram=False, legend_gap=15,
                       cmap='Blues', rasterized=True, show_colnames=False, show_rownames=True, row_names_side='right',
                       xlabel="", legend_hpad=0.2, xlabel_kws={'horizontalalignment':'center'},
                       xticklabels_kws={'labelrotation': 0, 'size': 14},
                        )

plt.xticks(rotation=0, ha='center')
plt.savefig(plot_out+"ComplexWoLin.pdf", bbox_inches='tight')
plt.close('All')


