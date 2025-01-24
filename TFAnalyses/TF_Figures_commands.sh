#! /bin/bash

# Commands to create TF-related figures and process tables.

## add info from extended TF table about target genes to then plot the oddsratio heatmap.
python getTargetGenes.py ../SNEEP_analysis/heatmapInput_oddsRatio_cutoffs_pseudocount_labeled.txt  ../SNEEP_analysis/TF_table_0810_extended.txt ../base_data/gencode.v38.annotation.gtf  heatmapInput_oddsRatio_cutoffs_pseudocount_labeled_targetGenes_11_10.txt
Rscript  ../SNEEP_analysis/heatmap_oddsRatio_2.0.R heatmapInput_oddsRatio_cutoffs_pseudocount_labeled_targetGenes_11_10.txt  oddsRatio_heatmap_11_10.pdf

## vulcano plot, piechart and barplot starnet data 
Rscript ../SNEEP_analysis/STARNET_DEG/vulcanPlot.R  ## results in a vulcanoPlot per tissue vulcanoPlot_TFs_*.pdf, piechart_allSTARNET_tissue.pdf, barplot_STARNET_tissue.pdf and a table of the sigTFs sigTFs_perTissue.txt

## figure 3D  -> similar for all coding and non-coding genes 
Rscript  heatmap_TF_SNV_REM.R ../figure_tables/2509_JointKnownGenes_All_Top20_TFSNVREM_Heatmap.tsv  heatmap_TF_SNV_REM_alleGenes_07_10.pdf
Rscript  heatmap_TF_SNV_REM.R ../figure_tables/2509_JointKnownGenes_non-codingRNA_Top20_TFSNVREM_Heatmap.tsv  heatmap_TF_SNV_REM_non_proteinCoding_07_10.pdf
Rscript  heatmap_TF_SNV_REM.R ../figure_tables/2509_JointKnownGenes_proteincoding_Top20_TFSNVREM_Heatmap.tsv  heatmap_TF_SNV_REM_proteinCodingGenes_07_10.pdf

## get number of protein coding and non-coding genes
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_ZNF610.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_ZNF610.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_ZFP14.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_ZFP14.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_ZIC4.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_ZIC4.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_ZBTB33.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_ZBTB33.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_TCFL5.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_TCFL5.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_MBD2.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_MBD2.txt
python getNumberncGenes_Fig3E.py  ../SNEEP_analysis/eQTL_TF_analysis/expressionPerTF/violinPlot_logFC_TFAP2E.txt ../base_data/gencode.v38.annotation.gtf  numberGenesBiotype_Fig3E_TFAP2E.txt

python getInfo_TFfunctional_enrichment.py ../figure_tables/ ../figure_tables/TFs_functionalEnrichment_detailed.txt  ## add manually the enriched pathway based on the name of orginal pathway -> ../figure_tables/TFs_functionalEnrichment_detailedPathways.txt
# get the detailedPathways table from excel file
## script that computes average p-value and fraction
python computeAverageFraction.py ../figure_tables/TFs_functionalEnrichment_detailedPathways.txt  ../figure_tables/TFs_functionalEnrichment_detailedPathways_average.txt

##plot again with adapted barPlot_TF_enrichment.R 
Rscript plotTFseQTLsEvidence.R  ../figure_tables/TFs_high_eQTL_evidence.txt ../figure_tables/TFs_functionalEnrichment_detailedPathways_average.txt  TFs_high_eQTL_evidence_25_06.pdf
