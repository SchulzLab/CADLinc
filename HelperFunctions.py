import pandas as pd
from pybedtools import BedTool
import gzip
from itertools import chain
import requests
import time
import math
from gprofiler import GProfiler
from matplotlib import cm
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import numpy as np
import PlottingFunctions

"""Collection of functions used for data processing."""


def gene_biotypes(gtf_file, gene_set=set()):
    """
    Gives a dict of Ensembl IDs with their biotype {g: {gtf: as noted in the gtf, general: manually summarized}.
    For the manually summarized biotypes, all gtf-biotypes containing 'RNA' or equal 'ribozyme' are converted to
    'non-coding RNA'. All biotypes containing 'pseudogene' are assigned 'pseudogene'. 'TcR gene' and 'Ig gene' are
    assigned to 'protein-coding'. Removes Ensembl ID version suffixes.

    Args:
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
    """

    def biotyper(biotype):
        """Translates the specific gtf biotypes into more general ones."""
        if 'pseudogene' in biotype:
            return 'pseudogene'
        elif 'RNA' in biotype or biotype == 'ribozyme':
            return "non-coding RNA"
        elif 'TR_' in biotype:
            return 'protein-coding'  #"TcR gene"  # T-cell receptor gene
        elif 'IG_' in biotype:
            return 'protein-coding'  #'Ig gene'  # Immunoglobulin variable chain gene
        elif biotype == 'protein_coding':
            return 'protein-coding'
        else:
            return biotype

    if gtf_file.endswith('.gz'):
        gtf_opener = gzip.open(gtf_file, 'rt')
    else:
        gtf_opener = open(gtf_file)
    biotype_dict = {}
    with gtf_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == 'gene':
                ensembl_id = entry.split('\t')[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                gene_name = entry[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                if not gene_set or ensembl_id in gene_set or gene_name in gene_set:
                    gtf_type = entry.split('\t')[8].split('gene_type "')[-1].split('"; ')[0]
                    biotype_dict[ensembl_id] = {'gtf': gtf_type,
                                                'general': biotyper(gtf_type)}
    return biotype_dict


def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS, also containing the number of transcripts, gene name and strand.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file.
        dict_only: Returns a dictionary instead of a BedTool's object.
        merge: If True, merges all overlapping promoters of the same gene into one row in the BedTool's object.
        open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
            BedTool's object. Can therefore be used to find genes whose promoter overlap a set of peaks, for example to
            find genes that are accessible.
    """
    if tss_type == '5':
        identifier = 'gene'
    elif tss_type == 'all':
        identifier = 'transcript'
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    tss_locs = {}
    with file_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == identifier:
                line = entry.strip().split('\t')
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    this_gene = line[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('";')[0]
                    if not gene_set or this_gene in gene_set or gene_name in gene_set:
                        if this_gene not in tss_locs:
                            tss_locs[this_gene] = {'chr': None, 'tss': set(), '#transcripts': 0}

                        tss_locs[this_gene]['chr'] = line[0]
                        tss_locs[this_gene]['name'] = gene_name
                        if line[6] == '+':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] > int(line[3])):
                                tss_locs[this_gene]['tss'] = {int(line[3])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[3]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '+'
                        if line[6] == '-':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] < int(line[4])):
                                tss_locs[this_gene]['tss'] = {int(line[4])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[4]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '-'

    if dict_only:
        return tss_locs

    promoter_bed = BedTool('\n'.join(chain(*[[vals['chr'] + '\t' + str(max([0, tss - int(extend) - 1])) + '\t' +
                                             str(tss + int(extend)) + '\t' + g + '\t.\t' + vals['strand'] for tss in vals['tss']]
                                             for g, vals in tss_locs.items()])), from_string=True)

    if open_regions and str(open_regions).lower() != "false":
        if type(open_regions) == str:
            open_regions = BedTool('\n'.join(['\t'.join(x.strip().split('\t')[:3]) for x
                                              in open(open_regions).readlines() if not x.startswith('#')]), from_string=True)
        promoter_bed = promoter_bed.intersect(open_regions)

    if merge:  # Flip the chr and Ensembl ID column to merge promoter of the same gene, and afterwards flip again.
        # Due to using the Ensembl ID as the first coordinate we don't need to consider the strand for merging (relying
        # on the annotations putting a gene only on one strand exclusively).
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True).sort().merge(c=[4, 5, 6], o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True)

    return promoter_bed


def gene_body_bed(gtf_file, gene_set=set(), dict_only=False):
    """
    From a gtf-file fetches the gene bodies, meaning start-end for the entries labelled as 'gene'.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        dict_only: Returns a dict {Ensembl ID: [chr, start, end, Ensembl ID, '.', strand]. Else, returns a bedtool-
            object with the coordinates per gene.
    """
    if gene_set:
        gene_set = set([g.split('.')[0] for g in gene_set])

    hits = []
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)

    with file_opener as gtf_in:
        for line in gtf_in:
            if not line.startswith('#') and line.split('\t')[2] == 'gene':
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line.strip().split('\t')[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    line = line.strip().split('\t')
                    gene = line[8].split('gene_id "')[-1].split('"; ')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('"; ')[0].split('.')[0]
                    if not gene_set or gene in gene_set or gene_name in gene_set:
                        hits.append([line[0], line[3], line[4], gene, line[5], line[6]])

    if dict_only:
        return {x[3]: x for x in hits}
    else:
        return BedTool('\n'.join(['\t'.join(x) for x in hits]), from_string=True)


def gene_location_bpwise(bed_dict, gtf_file, plot_path, tss_type='5', external_bed={}, palette='tab20',
                         formats=['pdf']):
    """
    Based on a gtf file builds bed-objects for Promoter (±200bp), Exons, UTRs and Introns, and then counts how many
    of the bp in the bed file(s) are located within those annotations and which are intergenic. All gene features are
    exclusive, overlaps are removed. Introns are gene bodies subtracted by all other features. The bed_dict can also
    be a list of bed files, to omit recreating the gtf-annotations each time. Creates a pie chart with
    the percentages in the given path. Also returns a dictionary with the bp-location of each bed-region and the
    total overlaps.

    Args:
        bed_dict: A dictionary with {title tag: bed-file or BedTools object}
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        tss_type: What to consider as promoter of genes, either '5' to use only the most 5' TSS, or 'all' to consider
            ±200bp around all annotated TSS of a gene.
        external_bed: An additional dictionary of bed-file or BedTools object which will be added as category for the
            intersection. This will be considered as highest priority, meaning the regions in there are removed from the
            gene-related features, and a bp overlapping external_bed will not be counted anywhere else. Multiple external
            bed-regions shouldn't overlap, that causes undefined outcomes.

    Returns:
        tuple:
            - **regions_locs**: For each entry in the bed_dict a dictionary with the bp-wise locations for each individual region.
            - **total_locs**: For each entry in the bed_dict the overall overlap of base pairs with each genomic feature.
    """
    if gtf_file.endswith('.gz'):
        bed_annotations = {'Promoter': gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       gzip.open(gtf_file, 'rt').readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     gzip.open(gtf_file, 'rt').readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    else:
        bed_annotations = {'Promoter': gene_window_bed(gtf_file, 200, tss_type=tss_type).sort().merge(),
                           'Exons': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                       open(gtf_file).readlines() if
                                                       not x.startswith('#') and x.split('\t')[2] == 'exon']),
                                            from_string=True).sort().merge(),
                           'UTR': BedTool('\n'.join(['\t'.join([x.strip().split('\t')[c] for c in [0, 3, 4]]) for x in
                                                     open(gtf_file).readlines() if
                                                     not x.startswith('#') and x.split('\t')[2] == 'UTR']),
                                          from_string=True).sort().merge()}
    # Remove the Promoter regions and UTRs from Exons and Promoter from UTRs, want to have exclusive annotations.
    bed_annotations['Exons'] = bed_annotations['Exons'].subtract(bed_annotations['Promoter']).subtract(
        bed_annotations['UTR'])
    bed_annotations['UTR'] = bed_annotations['UTR'].subtract(bed_annotations['Promoter'])

    introns_bed = BedTool(gtf_file).sort().merge()
    for remove in ['Promoter', 'Exons', 'UTR']:
        introns_bed = introns_bed.subtract(bed_annotations[remove])
    bed_annotations['Introns'] = introns_bed

    # Add potential external bed files and subtract it from all other regions.
    if external_bed:
        for external in external_bed:
            bed_annotations[external] = BedTool(external_bed[external]).sort().merge()
            for anno in bed_annotations:
                if anno not in external_bed:
                    bed_annotations[anno] = bed_annotations[anno].subtract(bed_annotations[external])

    regions_locs = {}
    total_locs = {}
    for tag, this_bed in bed_dict.items():
        if len(this_bed) == 0:
            print("Empty bed", tag)
            continue
        # Merge the bedfile to have unique bp, but keep track of what got merged and assemble it back later.
        bed = BedTool(this_bed).sort().merge(c=[1, 2, 3], o=['collapse'] * 3)
        org_beds = {'\t'.join(x.fields[:3]): [] for x in bed}
        for site in bed:
            org_beds['\t'.join(site.fields[:3])] += [
                '\t'.join([site.fields[3].split(',')[i], site[4].split(',')[i], site[5].split(',')[i]]) for i in
                range(site.fields[4].count(',') + 1)]
        bed_inter = {'\t'.join(x.fields[:3]): {a: 0 for a in bed_annotations.keys()} for x in bed}
        # Very broad marks rarely intersect only one annotation, so we count each bp.
        overall_inter = {a: 0 for a in list(bed_annotations.keys()) + ['Intergenic']}

        for annot, annot_bed in bed_annotations.items():
            intersection = bed.intersect(annot_bed, wo=True)
            for inter in intersection:
                bed_inter['\t'.join(inter.fields[:3])][annot] += int(inter.fields[-1])
                overall_inter[annot] += int(inter.fields[-1])

        # Intergenic is every base that was not assigned to any of the other annotations yet.
        for entry in bed_inter:
            missing_bp = abs(int(entry.split('\t')[2]) - int(entry.split('\t')[1])) - sum(bed_inter[entry].values())
            if missing_bp < 0:
                print("WARNING: There's a negative number of missing bp!")
            overall_inter['Intergenic'] += missing_bp
            bed_inter[entry]['Intergenic'] = missing_bp

        locations = {k: 0 for k in overall_inter.keys()}
        for annot in set(overall_inter.keys()):
            max_inter = len([x for x in bed_inter if max(bed_inter[x], key=bed_inter[x].get) == annot])
            print(annot, round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2), 'max', max_inter,
                  round(max_inter / len(bed) * 100, 2))
            locations[annot] = round(overall_inter[annot] / sum(overall_inter.values()) * 100, 2)

        loc_df = pd.DataFrame([[annot, locations[annot]] for annot in overall_inter.keys()],
                              columns=['Location', 'Overlap']).set_index('Location')
        PlottingFunctions.basic_pie(plot_df=loc_df, title=tag + '\n#' + str(len(BedTool(this_bed))), palette=palette,
                               numerate=False, legend_perc=True, formats=formats, legend_title='',
                               output_path=plot_path + (tag + "_GeneFeatureLocation_bpwise").replace(" ", '_'))

        # Now map the potentially merged regions back to all its original regions.
        org_inters = {}
        for inter, locs in bed_inter.items():
            for sub_i in org_beds[inter]:
                org_inters[sub_i] = locs
        regions_locs[tag] = org_inters
        total_locs[tag] = locations

    return regions_locs, total_locs



def go_enrichment(go_genes, title_tag='', out_tag='', max_terms='all', organism='hsapiens', background=None,
                  numerate=False, godf_only=False, wanted_sources=['GO:MF', 'GO:BP', 'KEGG', 'REAC', 'HP', 'WP'],
                  keywords={}, cmap='plasma', fig_width=None, fig_height=None, legend_out=None, rotation=45, font_s=16,
                  custom_dfs=None, formats=['pdf']):
    """Requires a dictionary with sets of genes. Will run GO enrichment for each of the sets. One plot will be written
    for each GO source. For only one gene set uses the x-axis for indicating the FDR, multiple sets will be
    separated on the x-axis and the FDR value shown as colour. Note, root terms of the databases are manually
    filtered out (e.g., HP root), but some might be missed. Uses the Python package of g:Profiler:
     https://biit.cs.ut.ee/gprofiler/gost.

    Args:
        go_genes: {class: [class genes] for class in classes}, or as sets. Safest identifiers are Ensembl IDs.
        max_terms: Allow a maximum of max_terms per dict key. Use 'all' to get all.
        organism: 'hsapiens' for human, 'mmusculus' for mouse. Check the g:Profiler page for all options.
        background: Same as go_genes but holding the set of background genes for each entry in go_genes. The keys of
            the two dictionaries will be matched. Leave empty to not use any background.
        numerate: Show the size of the gene sets in parentheses on the x-axis.
        godf_only: Skip the plotting step and only return the df.
        wanted_sources: Which databases to plot.
        keywords: A dictionary of {source: list of keywords}, e.g. {GO:BP: ['vascular', 'angio']}, to limit the
            plot to terms containing any of the listed strings. Capitalization doesn't matter, string comparison is done
            on the lowered strings.
        custom_dfs: Very specific use case. {gene set: {go_source: DF}}; In case an external df was used or
            a g:Profiler was customized, to just use the plotting part. Must have the same format as the g:Profiler DFs.

    Returns:
        - Returns a dict of {key df} with a df for each gene set as provided by the gprofiler package, which includes which genes matched to which terms.
    """
    keywords = {k: [t.lower() for t in vals] for k, vals in keywords.items()}  # We compare lower strings.
    cmap = cm.get_cmap(cmap)
    norm = plt.Normalize(0, 0.05)
    dict_keys = list(go_genes.keys())
    term_fetcher = {c: {s: [] for s in wanted_sources} for c in dict_keys}
    df_fetcher = {c: None for c in dict_keys}

    if not custom_dfs:
        gp = GProfiler(return_dataframe=True)
        for cell in dict_keys:
            if len(go_genes[cell]) >= 1:
                query_df = gp.profile(organism=organism, query=list(go_genes[cell]), no_evidences=False,
                                      background=None if not background else list(background[cell]))
                query_df = query_df[(~query_df['name'].str.contains("REACTOME")) & (query_df['name'] != 'WIKIPATHWAYS')
                                    & (~query_df['name'].str.contains("KEGG")) & (~query_df['name'].str.contains('HP root'))]
                query_df['Gene fraction'] = query_df['intersection_size'] / len(go_genes[cell])
                df_fetcher[cell] = query_df
                for source in wanted_sources:
                    term_fetcher[cell][source] = query_df[query_df['source'] == source]
            else:
                print(cell, 'empty gene set')
        if godf_only:
            return df_fetcher
    else:
        term_fetcher = custom_dfs

    for source in wanted_sources:
        print(source)
        term_collection = set()
        for cell in dict_keys:
            if len(term_fetcher[cell][source]) > 0:
                if str(max_terms).lower() == 'all':
                    term_collection = term_collection.union(set(term_fetcher[cell][source]['name'].to_list()))
                else:
                    term_collection = term_collection.union(set(term_fetcher[cell][source]['name'].to_list()[:max_terms]))

        these_terms = list(term_collection)
        if source in keywords:
            these_terms = [x for x in these_terms if np.any([k.lower() in x.lower() for k in keywords[source]])]
        # First count in how many cell types a term is present and sort according to that.
        term_occs = []
        for term in these_terms:
            hits = 0
            pvals = []
            for cell in dict_keys:
                curr_df = term_fetcher[cell][source]
                if len(curr_df) > 0:
                    if len(curr_df[curr_df['name'] == term]) == 1:
                        hits += 1
                        pvals.append(next(iter(curr_df[curr_df['name'] == term]['p_value'].values)))
            term_occs.append([term, hits, min(pvals)])
        sorted_terms = [x[0] for x in sorted(term_occs, key=lambda x: (x[1], -x[2]))]

        main_list = []  # X-pos, Y-pos, gene fraction, p-value.
        for x, cell in enumerate(dict_keys):
            num_genes = len(go_genes[cell])
            for y, term in enumerate(sorted_terms):
                curr_df = term_fetcher[cell][source]
                if len(curr_df) > 0:
                    match_entry = curr_df[curr_df['name'] == term]
                    if len(match_entry) == 1:
                        main_list.append([x, y, match_entry['intersection_size'].values[0] / num_genes,
                                          match_entry['p_value'].values[0]])
        if main_list:
            size_col = [s[2] for s in main_list]
            scaled_min, scaled_max = 20, 200
            if len(main_list) > 1 and min(size_col) != max(size_col):  # Otherwise all get the same size.
                scaled_sizes = [((scaled_max - scaled_min) * (s - min(size_col)) / (max(size_col) - min(size_col))) + scaled_min
                                for s in size_col]
            else:
                scaled_sizes = [scaled_max] * len(size_col)
            format_terms = []
            for term in sorted_terms:
                if len(term) > 40 and ' ' in term:
                    closest_blank = \
                    sorted([[i, abs(i - len(term) // 2)] for i, x in enumerate(term) if x == ' '], key=lambda x: x[1])[0][0]
                    format_terms.append(term[:closest_blank] + '\n' + term[closest_blank + 1:])
                elif len(term) > 40:
                    format_terms.append(term[:len(term) // 2] + '-\n' + term[len(term) // 2:])
                else:
                    format_terms.append(term)
            # If we have multiple gene groups, do a bubble plot aligning groups vertically.
            this_height = max([int(len(sorted_terms) * 0.4), 6]) if not fig_height else fig_height
            f, ax = plt.subplots(figsize=(fig_width if fig_width else (3 + len(dict_keys)), this_height))
            ax.yaxis.grid(True, which='major', color='#f0eded', zorder=1)
            ax.set_ylim(-0.5, len(these_terms) - 0.5)
            if len(go_genes) > 1:
                ax.set_xlim(-0.2, len(dict_keys) - 0.9)
                plt.scatter(x=[x[0] for x in main_list], y=[y[1] for y in main_list], c=[cmap(norm(c[3])) for c in main_list],
                            s=scaled_sizes, edgecolors=[cmap(norm(c[3])) for c in main_list], zorder=12)

                ax.set_xticks([x for x in range(len(dict_keys))])
                ax.set_xticklabels([c for c in dict_keys], size=font_s - 6, ha='right' if rotation else 'center',
                                       fontweight='bold', rotation=rotation)
                if numerate:
                    x_counts = {c: len(vals) for c, vals in go_genes.items()}
                    ax.set_xticklabels([x._text + '\n(#' + str(x_counts[x._text]) + ')' for x in ax.get_xmajorticklabels()])
                cax = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5)
                x_offset = cax.ax.get_position().x1
                cax.set_label('adj. p-value', size=font_s - 4)
                for c in range(len(dict_keys)):
                    plt.axvline(c, color='#c9c7c7', ymax=len(these_terms)-1, linewidth=0.5)
                ax.set_title(source + ": " + title_tag, size=font_s)

            else:  # If we only have one group, use the x-axis for the log-q value.
                ax.set_xlim(3, max([abs(np.log2(x[3])) for x in main_list])*1.1)
                plt.scatter(x=[abs(np.log2(x[3])) for x in main_list], y=[y[1] for y in main_list],
                            c=['#112791' for c in main_list],  # cmap(norm(c[3]))
                            s=scaled_sizes, edgecolors=['#112791' for c in main_list], zorder=12)
                ax.set_yticks([y for y in range(len(sorted_terms))])
                ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
                ax.set_xlabel('-log2 FDR', size=font_s - 4)
                ax.set_title(
                    source + ": " + title_tag + (' (#' + str(len(next(iter(go_genes.values())))) + ')') * numerate,
                    size=font_s)
                x_offset = 1.1

            ax.set_yticks([y for y in range(len(sorted_terms))])
            ax.set_yticklabels(format_terms, fontweight='bold', size=font_s - 8)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            if len(scaled_sizes) > 2:
                h1 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_min), color='black', linestyle='None')
                h2 = Line2D([0], [0], marker='o', markersize=np.sqrt((scaled_min + scaled_max) / 2), color='black',
                            linestyle='None')
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                leg = plt.legend([h1, h2, h3], [str(round(min(size_col), 4)), str(round((min(size_col) + max(size_col)) / 2, 4)),
                                          str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 6, title="Gene fraction", bbox_to_anchor=(legend_out if legend_out else x_offset+0.5, 1))
            else:
                h3 = Line2D([0], [0], marker='o', markersize=np.sqrt(scaled_max), color='black', linestyle='None')
                leg = plt.legend([h3], [str(round(max(size_col), 4))], loc="upper right", markerscale=1, scatterpoints=1,
                           fontsize=font_s - 6, title="Gene fraction", bbox_to_anchor=(legend_out if legend_out else x_offset+0.5, 1))
            leg.get_title().set_fontsize(11)
            if type(formats) != list:
                formats = [formats]
            for form in formats:
                f.savefig((out_tag + "_" + source + "_max"+str(max_terms)+"."+form).replace(' ', '').replace(':', ''),
                          bbox_inches='tight', format=form)
            plt.close()

        else:
            print("No enrichment", source)

    return df_fetcher


def match_genenames(gene_symbols, gtf_file, species='human', scopes="symbol, alias, uniprot", fields="ensembl, symbol"):
    """
    Takes a list of gene names to look for their matching Ensembl ID first in a gtf-file. Missing names will be
    queried via the API of https://mygene.info/. Only hits with a valid Ensembl ID will be returned. If multiple hits
    are found, the first one is used. Name mapping is fun.
    :param gene_symbols: List of gene names/symbols.
    :param: gtf_file: Gene annotation in gtf-style, can be gzipped.
    :param species: Give the name of the species, see below for the available ones.
    :param scopes: Where mygene.info will search for the gene symbols.
    :param fields: To which fields the output will be limited to.
    :return: mapped_names: Dict of {gene name: Ensembl ID} of the mappable names.
    :return: no_hits: Names which were not mappable via gtf-file nor mygene.info.
    """
    available_species = ['human', 'mouse', 'rat', 'fruitfly', 'nematode', 'zebrafish', 'thale-cress', 'frog' and 'pig']
    if species not in available_species:
        print("ERROR: species not available with name for mygene.info", species)
        print("Available are", available_species)
        return

    if gtf_file.endswith('.gz'):
        opener = gzip.open(gtf_file, 'rt')
    else:
        opener = open(gtf_file)
    name_id_map = {}
    for line in opener:
        if not line.startswith("#") and line.split('\t')[2] == 'gene':
            name_id_map[line.split('\t')[8].split('gene_name "')[-1].split('";')[0].lower()] = \
                line.split('\t')[8].split('gene_id "')[-1].split('";')[0].split('.')[0]

    gtf_missed = []
    mapped_names = {}
    for name in gene_symbols:
        if name.lower() in name_id_map:
            mapped_names[name] = name_id_map[name.lower()]
        else:
            gtf_missed.append(name)

    # Use the API from mygene for those we can't find in the gtf-file.
    if len(gtf_missed) > 0:
        print("Querying mygene for names", len(gtf_missed))
        query_data = {'species': species,
                      'scopes': scopes,
                      'fields': fields,
                      'ensemblonly': 'true'}
        query_n = len(gtf_missed)
        query_time = time.time()
        if query_n <= 1000:
            query_data['q'] = ' '.join(gtf_missed)
            res = requests.post('https://mygene.info/v3/query', query_data)
            res_json = res.json()
        else:
            # If the query is too long, we will need to break it up into chunks of 1000 query genes (MyGene.info cap)
            if query_n % 1000 == 0:
                chunks = query_n / 1000
            else:
                chunks = (query_n / 1000) + 1
            query_chunks = []
            for i in range(math.ceil(chunks)):
                start_i, end_i = i*1000, (i+1)*1000
                query_chunks.append(' '.join(gtf_missed[start_i:end_i]))
            res_json = []
            for chunk in query_chunks:
                query_data['q'] = chunk
                res = requests.post('https://mygene.info/v3/query', query_data)
                res_json = res_json+list(res.json())
        print('Batch query complete:', round(time.time()-query_time, 2), 'seconds')

        for entry in res_json[:len(gtf_missed)]:  # There can be appended entries that are just plain strings.
            if 'ensembl' in entry:
                if type(entry['ensembl']) == list:
                    mapped_names[entry['query']] = entry['ensembl'][0]['gene']  # Only keep the first hit.
                else:
                    mapped_names[entry['query']] = entry['ensembl']['gene']

    total_missed = [g for g in gene_symbols if g not in mapped_names]
    print("Non-matchable names", len(total_missed))
    return mapped_names, total_missed


def open_genes_multibed(bed_dict, annotation):
    """For each of the bed files in bed_dict {tag: bed file} returns the set of genes that intersect with at least
    one of their promoters with a peak."""
    promoter_bed = gene_window_bed(annotation, extend=200, tss_type='all')
    open_genes = {}
    for tag, bed in bed_dict.items():
        print(bed)
        if type(bed) == str:
            bed = BedTool(bed)
        open_genes[tag] = set([x.fields[3] for x in promoter_bed.intersect(bed, u=True)])

    return open_genes

