import os
import subprocess
import gzip
from pybedtools import BedTool
from liftover import get_lifter

"""
Takes peak files in one folder and merges potential replicates.
Based on the ABC-scoring creates a peak file which also lists the predicted target genes.
Enhancers without target genes are still included.

 Note: The whole processing is dependent on a fixed naming schema of the files:
 ID_Type_Column_GenomeVersion_ReplicateNumber.bed.
 Processing is done once for each ID, and needs to be referencable in the metadata file.
"""

metadata_file = "Enhancer_metadata.txt"
peak_folder = "peak_files/"
abc_folder = "ABC_scoring_avghg38_unique/"
inter_folder = "enhancer_interactions_avghg38_unique/"
if not os.path.isdir(abc_folder):
    os.mkdir(abc_folder)
if not os.path.isdir(inter_folder):
    os.mkdir(inter_folder)

# Input for the abc scoring and the lifting afterwards.
executable = "STARE/Dev_Code/STARE_ABCpp"
gene_annotation = "gencode.v38.annotation.gtf"
blocklisted_regions = "hg38-blacklist.v2.bed"
contact_folder = "ENCFF134PUN_avgHiC_hg38/"
# Read metadata to get the path to the Hi-C folder for each ID.
meta_head = {x: i for i, x in enumerate(open(metadata_file).readlines()[0].strip().split('\t'))}

# First, identify the peak files, and group them by ID.
already_scored = set([x.split('_')[0] for x in os.listdir(abc_folder)])
all_IDs = set([x.split('_')[0] for x in os.listdir(peak_folder) if
               not x.startswith('.') and len(x.split('_')) >= 4 and 'Duplicated' not in x]) - already_scored


# --------------------------------------------------------------------------------------------------------------
# Merging, Lifting, ABC-scoring
# --------------------------------------------------------------------------------------------------------------

def genome_lifter(region_list, converter):
    lifted_regions = []
    for region in region_list:
        chro, start, end = region[:3]
        hg38_start = converter[chro][int(start)]
        hg38_end = converter[chro][int(end)]
        # The conversion failed if one position is not liftable and when the strands differ. We add as condition
        # that the lifted region does not exceed a size of twice the original region.
        if len(hg38_start) > 0 and len(hg38_end) > 0 and hg38_start[0][2] == hg38_end[0][2] and abs(
                hg38_end[0][1] - hg38_start[0][1]) <= 2 * (int(end) - int(start)):
            if hg38_end[0][1] > hg38_start[0][1]:  # The conversion sometimes maps to the minus strand, mixing locs.
                lifted_regions.append([chro, str(hg38_start[0][1]), str(hg38_end[0][1])] + region[3:])
            else:
                # On the minus strand we have to +1 and reverse start and end to match the UCSC lift.
                lifted_regions.append([chro, str(hg38_end[0][1] + 1), str(hg38_start[0][1] + 1)] + region[3:])
    return lifted_regions


for cell in all_IDs:
    print(cell)
    matching_files = [x for x in os.listdir(peak_folder) if len(x.split('_')) >= 4 and x.split('_')[0] == cell]
    format_match = [x.split('_')[3].split('.')[0] == 'hg38' and len(x.split('_')) == 4 for x in matching_files]

    # Check if we already have a peak file in hg38 without replicates.
    if True not in format_match and len(matching_files) >= 1:
        print(matching_files)
        if len(format_match) > 1:  # Merge if we had multiple files.
            print('merging')
            merge_str = ''
            for rep in matching_files:
                merge_str += open(peak_folder + '/' + rep).read() + '\n'
            merged_bed = BedTool(merge_str, from_string=True).sort().merge(c=int(rep.split('_')[2].replace('c', '')),
                                                                           o='mean')
            # Add an ID column after the averaged activity.
            id_merged = '\n'.join([x + '\t' + str(i) for i, x in enumerate(str(merged_bed).strip().split('\n'))])
            out_file = '_'.join(rep.split('_')[:2]) + '_c4_' + rep.split('_')[3] + '.bed'
            open(peak_folder + '/' + out_file, 'w').write(
                id_merged)  # Even if it's hg38, one merged reference is useful.
        elif len(format_match) == 1:
            out_file = matching_files[0]

        if out_file.split('_')[3].split('.')[0] == 'hg38':
            final_file = out_file
        else:  # Lift the merged file or the hg38 to hg19.
            print('lifting')
            hg38_regions = genome_lifter([x.strip().split('\t') for x in open(peak_folder+'/'+out_file).readlines() if not x.startswith('#')],
                                         converter=get_lifter('hg19', 'hg38'))
            final_file = '_'.join(out_file.split('_')[:3]) + '_hg38.bed'
            open(peak_folder + '/' + final_file, 'w').write('\n'.join(['\t'.join(x) for x in hg38_regions]))

    elif True in format_match and len(matching_files) >= 1:
        final_file = matching_files[format_match.index(True)]  # We assume to only have one matching file.
    else:
        print('WARNING: Format not recognized', cell, matching_files)

    print(final_file)

    # Now we can call the ABC-soring on the potentially merged hg19 file.
    subprocess.call(
        executable + ' -b ' + peak_folder + '/' + final_file + ' -n ' + final_file.split('_')[2].replace('c', '') +
        ' -a ' + gene_annotation + ' -w 5000000 -f ' + contact_folder + ' -k 5000 -t 0.02 -o '
        + abc_folder + '/' + final_file.split('.')[0] + ' -c 12 -x ' + blocklisted_regions, shell=True)

    print('\n')


# --------------------------------------------------------------------------------------------------------------
# Writing the peak-interaction files.
# --------------------------------------------------------------------------------------------------------------

already_processed = set([x.split('_')[0] for x in os.listdir(inter_folder)])
abc_ids = set([x.split('_')[0] for x in os.listdir(abc_folder) if not x.startswith('.') and len(x.split('_')) >= 4]) - already_processed

for cell in abc_ids:
    # We need the scoredInteractions file as well as the original peak file to lift it. Both must be there.
    abc_file = [x for x in os.listdir(abc_folder) if x.split('_')[0] == cell and 'ABCpp_scoredInteractions' in x
                and x.split('_')[3].split('.')[0] == 'hg38' and x.endswith('.gz')][0]
    peak_file = [x for x in os.listdir(peak_folder) if x.split('_')[0] == cell and len(x.split('_')) == 4
                 and x.split('_')[3].split('.')[0] == 'hg38'][0]
    print(cell, abc_file, peak_file)

    # Now we map the genes onto enhancers.
    hg38_peaks_genes = {x.split('\t')[0].replace('chr', '') + ':' + x.split('\t')[1] + '-' + x.strip().split('\t')[2]: []
                        for x in open(peak_folder+'/'+peak_file).readlines() if not x.startswith('#')}
    with gzip.open(abc_folder+'/'+abc_file, 'rt') as abc_in:
        inter_head = {x: i for i, x in enumerate(abc_in.readline().strip().split('\t'))}
        for inter in abc_in:
            peak = inter.strip().split('\t')[inter_head['PeakID']]
            hg38_peaks_genes[peak].append([inter.strip().split('\t')[inter_head['Ensembl ID']].split('.')[0],
                                            float(inter.strip().split('\t')[inter_head['TSS-dist']])])

    hg38_tagged = [['chr' + x.split(':')[0], x.split(':')[-1].split('-')[0], x.split('-')[-1], x] for x in hg38_peaks_genes]

    enhancer_genes = []
    missing_genes = set()
    for peak in hg38_tagged:
        gene_hits = []
        for gene in hg38_peaks_genes[peak[3]]:
            gene_hits.append(gene[0])
        if gene_hits:
            gene_str = ','.join(gene_hits)
        else:
            gene_str = '-'
        enhancer_genes.append(peak + [gene_str])

    enhancer_gene_bed = BedTool('\n'.join(['\t'.join(x) for x in enhancer_genes]), from_string=True)

    enhancer_out = inter_folder + '/' + '_'.join(peak_file.split('_')[:3]) + '_hg38_EnhancerInteractions.bed'
    print(enhancer_out, '\n')
    open(enhancer_out, 'w').write(str(enhancer_gene_bed))

