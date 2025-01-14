import pandas as pd
from itertools import chain
import os
import numpy as np
import gzip
import pickle
import itertools
from collections import Counter
from pybedtools import BedTool
import HelperFunctions
import PlottingFunctions


"""Collect the multitude of data to define the gene set we get from SNEEP after setting specific filters, and
join that with the GATES genes from an external file (see GATES_FDR.py), to produce a variety of plots."""

tag = '2509_JointKnownGenes'
plot_out = "Plots/" + tag + '_'
figure_tables_out = "FigureTables/" + tag + "_"
godf_out = "GO_DFs/" + tag + "_"

meta_file = "Enhancer_metadata.txt"
annotation = "gencode.v38.annotation.gtf.gz"
gene_name_map = {x.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]:
                     x.split('\t')[8].split('gene_name "')[-1].split('"; ')[0]
                 for x in gzip.open(annotation, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'gene'}
gene_id_map = {x.split('\t')[8].split('gene_name "')[-1].split('"; ')[0]:
                     x.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                 for x in gzip.open(annotation, 'rt').readlines() if not x.startswith('#') and x.split('\t')[2] == 'gene'}
gene_biotype_map = HelperFunctions.gene_biotypes(annotation)
gene_lengths = HelperFunctions.gene_body_bed(annotation, dict_only=True)
gene_lengths = {g: int(vals[2]) - int(vals[1]) for g, vals in gene_lengths.items()}
gene_bodies = HelperFunctions.gene_body_bed(gtf_file=annotation, dict_only=True)

sneep_file = 'sneep_cadlincGWAS_2.52x0.00001_22_05.txt'
ld_file = "summary_filtered_betaSwapped_2.52x0.00001_hg38_22_02_lead_proxySNPs_swapped_uniq_allLDSNVs.txt"
gene_based_file = "result_gene_based_test_genebody_Filtered.txt"
excluded_genes_file = "ExcludedGenes.txt"
excluded_genes = set([x.strip().split('\t')[0] for x in open(excluded_genes_file).readlines()])
nc_conserved_file = "rbbh_gene_alignments.csv"
biomart_orth_file = "BioMart_GRCh38.p14_MouseOrthologues.txt.gz"
coexpression_file = 'FisherTestCoExpressedGenes_sneep_geneBased_100_greater_JointKnownGenes_2509.txt'
gwas_file = "FINAL_1MH_CAD_GWAS_summary_stats.tsv.gz"
vcf_out = "VEP"
all_tfsnv_file = "positionAllRegSNPs_sorted.bed.gz"
all_tfsnvs = set([x.strip().split('\t')[0].replace('chr', '') + ':' + x.strip().split('\t')[2] for x in gzip.open(all_tfsnv_file, 'rt').readlines()])

all_rems_file = "AllREMsMerged_200624.bed"
enhancer_folder = "enhancer_interactions_avghg38_unique/"

# Files of other gene sets.
gwas_gene_file = "250924_JointKnownGenes_EnsemblIDs.txt"
colocstar_gene_file = "starnet.alltss.Coloc_Gene_.abf.res.pp4_0.5.xlsx"
colocgtex_file = "gtex.alltss.Coloc_Gene_.abf.res.pp4_0.5.xlsx"
cad_loci_genes = set([x.strip() for x in open(gwas_gene_file).readlines()]) - excluded_genes
coloc_gtex_genes = set(pd.read_excel(colocgtex_file)['geneid']) - excluded_genes
coloc_starnet_genes = set(pd.read_excel(colocstar_gene_file)['gene']) - excluded_genes
nc_conserved = set([x.split('.')[0] for x in pd.read_table(nc_conserved_file, sep='\t', header=0)['query_gene_y'].values])

# Get the Schnitzler genes.
schnitzler_file = "41586_2024_7022_MOESM4_ESM.xlsx"
schnitzler_studies = pd.read_excel(schnitzler_file, header=0, engine='openpyxl', sheet_name='Suppl.Table.17')
schnitzler_studies = schnitzler_studies[[c for c in schnitzler_studies.columns if not c.startswith('V2G2P') and 'EC' not in c and c not in ['known_CAD_gene', 'Li_TWAS']]]
schnitzler_id_map, missed = HelperFunctions.match_genenames(schnitzler_studies['Gene'].values, annotation, species='human')
schnitzler_studies = schnitzler_studies[schnitzler_studies['Gene'].isin(schnitzler_id_map)]
schnitzler_studies['Ensembl ID'] = [schnitzler_id_map[g] for g in schnitzler_studies['Gene'].values]
schnitzler_sets = {c: set(schnitzler_studies['Ensembl ID'][map(bool, schnitzler_studies[c])]) - excluded_genes for c in schnitzler_studies.columns if c not in ['Ensembl ID', 'Gene']}

# Take the gene-based test genes.
gene_based_df = pd.read_table(gene_based_file, header=0, sep='\t')
gene_based_genes = set(gene_based_df[gene_based_df['Filtered'] == True]['gene_list...1.'].values)
gene_based_dict = gene_based_df.set_index('gene_list...1.').to_dict(orient='index')

# Get the SNPs in LD, which we can later map via the rsIDs. Also fetch the hg19 position and nonrisk/risk allele from that.
snp_lds = {}
snp_hg19_map = {}
with open(ld_file) as ld_in:
    ld_header = {x: i for i, x in enumerate(ld_in.readline().strip().split('\t'))}
    for row in ld_in:
        row = row.strip().split('\t')
        this_snp = row[ld_header['rsId']]
        if this_snp != '-' and this_snp != '.':
            snp_map_id = this_snp
            if this_snp not in snp_lds:
                snp_lds[this_snp] = set(row[ld_header['ldSNVs']].split(','))
            else:
                snp_lds[this_snp] = snp_lds[this_snp].union(set(row[ld_header['ldSNVs']].split(',')))
        else:  # Use the hg38 position as identifier in case of missing rsID or self-made XX ID.
            snp_map_id = row[ld_header['chr']] + ":" + row[ld_header['start']] + "-" + row[ld_header['end']]
        snp_hg19_map[snp_map_id] = {'chr_hg19': row[ld_header['position_hg19']].split(':')[0].replace('chr', ''), 'POS_hg19': row[ld_header['position_hg19']].split('-')[1],
                                    'Allele1': row[ld_header['allele1']], 'Allele2': row[ld_header['allele2']],
                                    'chr_hg38': row[ld_header['chr']].replace('chr', ''), 'POS_hg38': row[ld_header['end']]}

sneep_df = pd.read_table(sneep_file, header=0, sep="\t")

# Get an ordering of the cell types and a matching class for each.
rem_meta_df = pd.read_table(meta_file, sep='\t', header=0)
cell_types = rem_meta_df['ID'].values
cell_class = {x[0]: x[1] for x in rem_meta_df[['ID', 'Class']].values}
cell_labels = {x[0]: x[1] for x in rem_meta_df[['ID', 'Label']].values}
lineage_cells = {l: set() for l in cell_class.values()}
for cell in cell_class:
    lineage_cells[cell_class[cell]].add(cell)


peak_files = {x.split('_')[0]: enhancer_folder + '/' + x for x in os.listdir(enhancer_folder) if 'EnhancerInteractions' in x}
open_genes_pkl = "OpenGenesCell.pkl"
if not os.path.isfile(open_genes_pkl):
    open_genes = HelperFunctions.open_genes_multibed(peak_files, annotation)
    pickle.dump(open_genes, open(open_genes_pkl, 'wb'))
else:
    open_genes = pickle.load(open(open_genes_pkl, 'rb'))
all_open = set().union(*open_genes.values())

# ----------------------------------------------------------------------------------------------------------------
# Get affected enhancers and their ABC-score and do loads of filtering
# ----------------------------------------------------------------------------------------------------------------
# Get the set of enhancer strings that intersect a regSNP from SNEEP.
sneep_genes = set(chain(*[x.split(',') for x in sneep_df['ensemblIDs'].values])) - set('.') - excluded_genes
snp_rems = set(chain(*[x.split(',') for x in sneep_df['REM_positions'].values])) - {'.'}
gene_regrems = {g: {'regSNPs': set(), 'regREMs': set(), 'rsIDs': set(),
                    'cell_rsIDs': {c: set() for c in cell_types}} for g in sneep_genes}
cell_gene_regrems = {c: {g: set() for g in sneep_genes} for c in cell_types}  # {Cell: Gene} works better later.
snp_tf_map = {rs: {} for rs in snp_hg19_map}
snp_rem_map = {c: {} for c in cell_types}
snp_pos = {}
regrem_collector = ""  # For faster merging afterwards, we replace the chr with the Ensembl ID, we stay within a chr.
for entry in sneep_df.to_dict(orient='records'):
    for gene, rem, cell in zip(entry['ensemblIDs'].split(','), entry['REM_positions'].split(','), entry['REMIds'].split(',')):
        if gene in sneep_genes:
            if gene in open_genes[cell]:  # Only consider CT-specific regREMs if a gene was open.
                regrem_collector += gene + '\t' + rem.split(':')[1].split('-')[0] + '\t' + rem.split('-')[1] + '\n'
                cell_gene_regrems[cell][gene].add(rem)
            gene_regrems[gene]['regSNPs'].add(entry['SNP_position'])
            snp_pos[entry['rsID']] = entry['SNP_position']
            gene_regrems[gene]['cell_rsIDs'][cell].add(entry['rsID'])
            gene_regrems[gene]['rsIDs'].add(entry['rsID'])
            snp_rem_map[cell][entry['rsID']] = rem
            snp_tf_map[entry['rsID']][entry['TF']] = str(entry['log_pvalueBindAffVar1_pvalueBindAffVar2'])
regrem_counts = {c: {g: len(val) for g, val in c_vals.items()} for c, c_vals in cell_gene_regrems.items()}
# Now merge the regREMs across cell types.
merged_regrems = BedTool(regrem_collector, from_string=True).sort().merge()
for merg in merged_regrems:
    gene_regrems[merg.fields[0]]['regREMs'].add(merg.fields[1] + '-' + merg.fields[2])

# Find genes that have ≥2 regREMs with non-LD SNPs in a cell type where it's open.
genes_passed_filters = set()
filtered_genes_cell = {c: set() for c in cell_types}
for cell in regrem_counts:
    for gene in sneep_genes:
        if regrem_counts[cell][gene] >= 2 and gene in open_genes[cell]:
            rs_combos = list(itertools.combinations(gene_regrems[gene]['cell_rsIDs'][cell], 2))
            for combo in rs_combos:  # Check each pair of SNPs whether they are in LD.
                if snp_rem_map[cell][combo[0]] != snp_rem_map[cell][combo[1]]:
                    if combo[0] not in snp_lds or combo[1] not in snp_lds:
                        genes_passed_filters.add(gene)  # Means that we have no SNPs in LD for one of them.
                    else:
                        if combo[1] not in snp_lds[combo[0]] and combo[0] not in snp_lds[combo[1]]:
                            genes_passed_filters.add(gene)  # We only need one pair of SNPs that is not in LD.
                            filtered_genes_cell[cell].add(gene)

final_gene_set = genes_passed_filters.union(gene_based_genes)

# ----------------------------------------------------------------------------------------------------------------
# Specificity of gene expression from IHEC.
# ----------------------------------------------------------------------------------------------------------------
expression_matrix_file = "genes_TPM.csv.gz"
ihec_meta_file = "IHEC_metadata_harmonization.v1.1.csv"
experiment_meta_file = "230723_epiatlas_metadata.csv"
sample_meta = pd.read_csv(ihec_meta_file, header=0)
experiment_df = pd.read_csv(experiment_meta_file, header=0)
meta_df = pd.merge(experiment_df[['uuid', 'epirr_id_without_version']], sample_meta, on="epirr_id_without_version", how='inner')
uuid_cell_map = {x['uuid']: x['harmonized_sample_ontology_intermediate'] for x in meta_df.to_dict(orient='records')}

expression_df = pd.read_table(expression_matrix_file, sep=',', header=0, index_col='id_col')
expression_df.columns = [uuid_cell_map[c.split('.')[-1]] for c in expression_df.columns]
expression_df = expression_df.groupby(by=expression_df.columns, axis=1).mean()
ubi_expression = {k.split('.')[0]: val/expression_df.shape[1] for k, val in (expression_df >= 0.5).sum(axis=1).items()}

# ----------------------------------------------------------------------------------------------------------------
# Paper table
# ----------------------------------------------------------------------------------------------------------------
gene_table = pd.DataFrame([[g, gene_name_map[g], g in cad_loci_genes, g in coloc_gtex_genes, g in coloc_starnet_genes] for g in final_gene_set],
                               columns=['Ensembl ID', 'Gene name', 'known CAD GWAS gene', 'Coloc GTEx', 'Coloc STARNET'])
gene_table['found via'] = [', '.join(['epigenome'*(g in genes_passed_filters), 'GATES'*(g in gene_based_genes)]).strip(', ') for g in gene_table['Ensembl ID'].values]

# Add whether a gene was found in one of the Schnitzler studies.
for schnitz in schnitzler_sets:
    gene_table['Schnitzler: '+schnitz] = [g in schnitzler_sets[schnitz] for g in gene_table['Ensembl ID'].values]

# Add whether they are coexpressed with other known CAD loci genes.
coexpression = pd.read_table(coexpression_file, sep='\t', header=0).set_index("gene").to_dict(orient='index')
coexpress_labels = []
for gene in gene_table['Ensembl ID'].values:
     if gene not in coexpression:
         coexpress_labels.append("insufficient expression data")
     else:
         if coexpression[gene]['FDR_corrected'] == '-':
             coexpress_labels.append("insufficient co-expressed protein coding genes")
         elif float(coexpression[gene]['FDR_corrected']) <= 0.05:
             coexpress_labels.append("True")
         elif float(coexpression[gene]['FDR_corrected']) > 0.05:
             coexpress_labels.append("False")
gene_table['Co-expressed with CAD GWAS genes'] = coexpress_labels
gene_table['Expression breadth'] = [None if g not in ubi_expression else ubi_expression[g] for g in gene_table['Ensembl ID'].values]
gene_table['Biotype'] = [gene_biotype_map[g]['gtf'] for g in gene_table['Ensembl ID'].values]
gene_table['General biotype'] = [gene_biotype_map[g]['general'] for g in gene_table['Ensembl ID'].values]

for label, idx in [['Chr', 0], ['Start', 1], ['End', 2], ['Strand', 5]]:
    gene_table[label] = [gene_bodies[g][idx] for g in gene_table['Ensembl ID'].values]

genes_open_cells = {g: set() for g in all_open}
for cell in open_genes:
    for g in open_genes[cell]:
        genes_open_cells[g].add(cell)
gene_table['Cell types with open promoter'] = [','.join(genes_open_cells[g]) for g in gene_table['Ensembl ID'].values]

# Add whether a gene is conserved.
biomart_orth_df = pd.read_table(biomart_orth_file, sep='\t', header=0)
biomart_conserved = set(biomart_orth_df[(biomart_orth_df['Mouse orthology confidence [0 low, 1 high]'] > 0)]['Gene stable ID'])
conserved_biotypes = Counter(['not found' if g not in gene_biotype_map else gene_biotype_map[g]['general'] for g in biomart_conserved])
gene_table['Conserved in mouse'] = [g in biomart_conserved or g in nc_conserved for g in gene_table['Ensembl ID'].values]

gene_table['Cell types with open promoter and ≥2 CREs with non-LD TF-SNVs'] = [','.join([c for c in cell_types if g in filtered_genes_cell[c]]) for g in gene_table['Ensembl ID'].values]
gene_table['Number ≤1% FDR TF-SNVs'] = ['' if g not in gene_regrems else len(gene_regrems[g]['rsIDs']) for g in gene_table['Ensembl ID'].values]
gene_table['Number affected TFs'] = [0 if g not in gene_regrems else len(set().union(*[snp_tf_map[rs] for rs in gene_regrems[g]['rsIDs']])) for g in gene_table['Ensembl ID'].values]
gene_table['Number merged CREs with ≤1% FDR TF-SNV from cell types with open promoter'] = [0 if g not in gene_regrems else len(gene_regrems[g]['regREMs']) for g in gene_table['Ensembl ID'].values]
gene_table['GATES FDR'] = ["" if g not in gene_based_dict else gene_based_dict[g]['FDR'] for g in gene_table['Ensembl ID'].values]

gene_table.to_csv(tag+'_GeneTable.txt', sep='\t', header=True, index=False)

# ----------------------------------------------------------------------------------------------------------------
# merged regREM distributions
# ----------------------------------------------------------------------------------------------------------------
sneep_sub_df = gene_table[gene_table['Ensembl ID'].isin(genes_passed_filters)]
sneep_sub_df['Number of CREs with TF-SNVs'] = sneep_sub_df['Number merged CREs with ≤1% FDR TF-SNV from cell types with open promoter']
PlottingFunctions.basic_hist(sneep_sub_df, x_col='Number of CREs with TF-SNVs', element='bars',
                        bin_num=21, palette=['#FFC20A', '#0C7BDC'], stat='proportion', hue_col='General biotype', hue_order=['protein-coding', 'non-coding RNA'],
                        title="", output_path=plot_out+"subPCncRNATestFillBar", alpha=1, xsize=9, ysize=5,
                        legend_title=False, font_s=14, fill=True, edgecolour='black', multiple='dodge', shrink=0.85)

# ----------------------------------------------------------------------------------------------------------------
# Plot expression ubiquitousness
# ----------------------------------------------------------------------------------------------------------------
PlottingFunctions.basic_hist(gene_table, x_col='Expression breadth', title="", output_path=plot_out, hue_order=['protein-coding', 'non-coding RNA'],
                        xsize=9, ysize=5, palette=['#FFC20A', '#0C7BDC'], legend_out=0.6, legend_title=False, fill=True,
                        alpha=0.6, bin_num=30, stat='proportion', hue_col='General biotype', font_s=14, element='bars')

# ----------------------------------------------------------------------------------------------------------------
# Venn of the gene sets.
# ----------------------------------------------------------------------------------------------------------------
gene_set_collections = {"candidate CAD genes": set(final_gene_set),
                        'known CAD GWAS genes': cad_loci_genes,
                        'Coloc GTEx': coloc_gtex_genes,
                        'Coloc STARNET': coloc_starnet_genes}
PlottingFunctions.upset_plotter(gene_set_collections, max_groups=None, sort_categories_by='input', y_label='Shared genes', element_size=48,
                           title_tag='Shared genes', show_percent=False, plot_path=plot_out+"Genesets", font_enhancer=6, intersection_plot_elements=8)
PlottingFunctions.venn_from_list([set(final_gene_set), cad_loci_genes, coloc_gtex_genes.union(coloc_starnet_genes)], ['', '', ''], plot_path=plot_out+"VennGeneSets",
                   scaled=True, linestyle='', number_size=16, blob_colours=['#9A201E', '#EEB576', '#3475C4'])

# ----------------------------------------------------------------------------------------------------------------
# Intersection known CAD genes and biotypes
# ----------------------------------------------------------------------------------------------------------------
cad_loci_pcgenes = set([g for g in cad_loci_genes if g in gene_biotype_map and gene_biotype_map[g]['general'] == 'protein-coding'])
gene_biotypes = {'All genes': Counter([vals['general'] for vals in gene_biotype_map.values()])}
for set_tag, gene_set, g_background in [['candidate CAD genes', final_gene_set, all_open - excluded_genes],
                                        ['Coloc GTEx', coloc_gtex_genes, set(gene_name_map.keys())],
                                        ['Coloc STARNET', coloc_starnet_genes, set(gene_name_map.keys())]]:
    gene_biotypes[set_tag] = Counter([gene_biotype_map[g]['general'] for g in gene_set if g in gene_biotype_map])

gene_biotypes['known CAD GWAS genes'] = Counter([gene_biotype_map[g]['general'] for g in cad_loci_genes if g in gene_biotype_map])
gene_biotype_df = pd.DataFrame.from_dict(gene_biotypes, orient='index')
PlottingFunctions.stacked_bars(plot_df=gene_biotype_df.loc[['candidate CAD genes', 'known CAD GWAS genes', 'Coloc GTEx', 'Coloc STARNET']][['protein-coding', 'non-coding RNA']],
                          x_col='Gene set', y_cols=['protein-coding', 'non-coding RNA'], y_label='genes', palette=['#FFC20A', '#0C7BDC'],
                          title="", output_path=plot_out+"Biotypes", x_size=5.2, y_size=4, legend_out=1.65,
                          rotation=0, legend=True, fraction=True, numerate=True, vertical=True)

# Get the exact biotype numbers for candidates and uniquely novel genes.
candidates_biotypes = Counter(gene_table['Biotype'])
unique_novel = final_gene_set - cad_loci_genes - coloc_gtex_genes - coloc_starnet_genes
novel_biotypes = Counter(gene_table[gene_table['Ensembl ID'].isin(unique_novel)]['Biotype'])
nonnovel_biotypes = Counter(gene_table[~gene_table['Ensembl ID'].isin(unique_novel)]['Biotype'])
gtf_biotype_df = pd.DataFrame.from_dict({"candidate\nCAD genes": candidates_biotypes, 'candidate CAD genes\nwith evidence': nonnovel_biotypes, 'novel candidate\nCAD genes': novel_biotypes}).T
gtf_biotype_df.columns = [c.replace('_', ' ') for c in gtf_biotype_df.columns]
PlottingFunctions.stacked_bars(plot_df=gtf_biotype_df, x_col='Gene set', y_cols=gtf_biotype_df.columns, y_label='genes',
                          palette=ColoursAndShapes.glasbey_cool, title="", output_path=plot_out+"CandidatesGTFBiotypes",
                          x_size=9, y_size=6, legend_out=1.35, rotation=0, legend=True, fraction=True, numerate=True,
                          vertical=False)

# ----------------------------------------------------------------------------------------------------------------
# Tables of highest affected genes per cell type.
# ----------------------------------------------------------------------------------------------------------------
# Rank by total number of regREMs across cell types.
regrem_df = pd.DataFrame.from_dict(regrem_counts, orient='index')[genes_passed_filters]
regrem_df['Class'] = [cell_class[c] for c in regrem_df.index]
regrem_df['Cell type'] = [cell_labels[c] for c in regrem_df.index]
ranks = 20
for filter_tag, filter_set in [['All', genes_passed_filters], ['protein coding', [g for g in genes_passed_filters if gene_biotype_map[g]['general'] == 'protein-coding']],
                               ['non-coding RNA', [g for g in genes_passed_filters if gene_biotype_map[g]['general'] == 'non-coding RNA']]]:
    top_genes = list(sorted(filter_set, key=lambda x: (len(gene_regrems[x]['regREMs']), len(gene_regrems[x]['rsIDs']),
                                                       len([c for c in cell_types if x in filtered_genes_cell[c]])), reverse=True))[:ranks+1]
    aggregate_df = regrem_df[top_genes].T
    aggregate_df.columns = [cell_labels[c] for c in aggregate_df.columns]
    aggregate_df = gene_table.join(aggregate_df, on='Ensembl ID', how='right')
    if filter_tag == 'non-coding RNA':  # Merge two MIRs at the same locus into one name and drop the other.
        aggregate_df = aggregate_df[aggregate_df['Gene name'] != 'MIR212']
        aggregate_df.loc[aggregate_df['Gene name'] == 'MIR132', 'Gene name'] = 'MIR132/MIR212'
    # Reduce size again
    aggregate_df = aggregate_df.iloc[:ranks]
    aggregate_df.to_csv(figure_tables_out + filter_tag.replace(' ', '') + "_Top"+str(ranks)+"_TFSNVREM_Heatmap.tsv",
                        header=True, sep='\t', index=False)

# ----------------------------------------------------------------------------------------------------------------
# GO of affected genes
# ----------------------------------------------------------------------------------------------------------------
# First enrichment of all the filtered genes.
keep_sources = ['GO:MF', 'GO:BP', 'KEGG', 'REAC', 'HP', 'WP']
final_genes_go = HelperFunctions.go_enrichment({"identified CAD genes": final_gene_set}, title_tag="Candidate CAD genes", out_tag=plot_out + "JointDG",
                           max_terms='all', font_s=16, numerate=True, organism='hsapiens', cmap='plasma', background={"identified CAD genes": all_open - excluded_genes})
final_genes_go['identified CAD genes'].to_csv(godf_out + "CandidateCADGenes_gProfiler.tsv", sep='\t', header=True, index=False)

# Manually selected among the top 20.
wp_terms = ['Cholesterol metabolism', 'Canonical and non canonical TGF B signaling', 'Familial hyperlipidemia type 3',
            'Familial hyperlipidemia type 1', 'Lipid particles composition', 'Familial hyperlipidemia type 5',
            'Fatty Acids and Lipoproteins Transport in Hepatocytes', 'Familial hyperlipidemia type 2',
            'Familial hyperlipidemia type 4', 'TGF beta signaling pathway', 'VEGFA VEGFR2 signaling',
            'Interleukin 11 signaling pathway', 'Triacylglyceride synthesis', 'miRNA targets in ECM and membrane receptors']
intermed_go = HelperFunctions.go_enrichment({"identified CAD genes": final_gene_set}, title_tag="identified CAD genes", out_tag=plot_out + "JointDGFilterWP",
                               max_terms=20, font_s=16, numerate=True, organism='hsapiens', cmap='plasma', wanted_sources=['WP'], rotation=0,
                               background={"identified CAD genes": all_open - excluded_genes}, keywords={"WP": wp_terms})

# And even another version where we remove redundant terms, altering the df from the previous function call.
wp_df = intermed_go['identified CAD genes']
wp_df = wp_df[(wp_df['source'] == 'WP') & (wp_df['name'].isin(wp_terms)) & (wp_df['name'] != 'TGF beta signaling pathway')]
fam_hyperlip = wp_df[wp_df['name'].str.contains('Familial hyperlipidemia')]
avg_fam_hyperlip = fam_hyperlip.mean(axis=0, numeric_only=True)
for col in fam_hyperlip.columns:
    if fam_hyperlip.dtypes[col] != 'float64':
        avg_fam_hyperlip[col] = fam_hyperlip.iloc[0][col]
avg_fam_hyperlip['name'] = 'Familial hyperlipidemia'
wp_df = wp_df[~wp_df['name'].str.contains('Familial hyperlipidemia')]
wp_df = pd.concat([wp_df, pd.DataFrame(avg_fam_hyperlip).T], axis=0)
_ = HelperFunctions.go_enrichment({"identified CAD genes": final_gene_set}, title_tag="identified CAD genes", out_tag=plot_out + "JointDGCustomWP",
                               max_terms=20, fig_width=2.5, fig_height=3, font_s=16, numerate=True, organism='hsapiens', cmap='plasma', wanted_sources=['WP'], rotation=0,
                               background={}, custom_dfs={"identified CAD genes": {"WP": wp_df}})

# ----------------------------------------------------------------------------------------------------------------
# Per lineage number of affected genes.
# ----------------------------------------------------------------------------------------------------------------
affected_cell = filtered_genes_cell
lineages = list(dict.fromkeys(cell_class.values()))  # To maintain the order.
lineage_genes = pd.DataFrame([[cell_labels[c], len(filtered_genes_cell[c]), cell_class[c]] for c in cell_types],
                                columns=['Cell type', 'Number of open genes with\n≥ 2 CREs with non-LD TF-SNVs', 'Lineage'])
lineage_means = {x: np.mean(lineage_genes[lineage_genes['Lineage'] == x]['Number of open genes with\n≥ 2 CREs with non-LD TF-SNVs'].values) for x in lineages}

# And split by biotype, for which we need a separate DF.
lineage_genes_biotypes = pd.DataFrame([[cell_labels[c], len([g for g in filtered_genes_cell[c] if gene_biotype_map[g]['general'] == 'protein-coding']), cell_class[c], 'protein-coding'] for c in cell_types] +
                                [[cell_labels[c], len([g for g in filtered_genes_cell[c] if gene_biotype_map[g]['general'] == 'non-coding RNA']), cell_class[c], 'non-coding RNA'] for c in cell_types],
                                columns=['Cell type', 'Number of genes\nregulated by TF-SNVs', 'Lineage', 'General biotype'])
PlottingFunctions.basic_violin(lineage_genes_biotypes, y_col='Number of genes\nregulated by TF-SNVs',
                         x_col='Lineage', x_order=[k for k, v in sorted(lineage_means.items(), key=lambda x: x[1], reverse=True)],
                         title="", hue_col='General biotype', output_path=plot_out, numerate=True, jitter_colour=['#FFC20A', '#0C7BDC'],
                         xsize=16, ysize=7, rotation=90, font_s=14, jitter=True, boxplot=True, boxplot_meanonly=False, saturation=1,
                         palette=['#FFC20A', '#0C7BDC'], legend_title=False, jitter_size=8, vertical_grid=True)



