# CADLinc

Repository to document scripts and software used in the manuscript "Cell type-specific epigenetic regulatory 
circuitry of coronary artery disease loci". This repository is not meant to be ready-to-use software to be applied to
new data, but as documentation of the manuscript analyses.


- [DataIntegration](https://github.com/SchulzLab/CADLinc/tree/main/DataIntegration/): Scripts that combine the information of other analyses, including plot creation for Figures 1 and 2.
- [SNEEPAnalyses](https://github.com/SchulzLab/CADLinc/tree/main/SNEEPAnalyses/): Gathering of SNV information, background sampling and running the SNEEP pipeline. Follow the commands in workflow_3.0.sh.
- [TFAnalyses](https://github.com/SchulzLab/CADLinc/tree/main/TFAnalyses/): Analyses and plotting on TF-level (Figure 3). The commands are written in TF_Figures_command.sh.
- [GATES](https://github.com/SchulzLab/CADLinc/tree/main/GATES/): How the GATES test was run, follow the commands in commands.sh.
- [PheWAS_Coloc](https://github.com/SchulzLab/CADLinc/tree/main/PheWAS_Coloc/): GTEx and STARNET data processing and PheWAS, including figure generation (Figure 4).
- [GenomicTrackView](https://github.com/SchulzLab/CADLinc/tree/main/GenomicTrackView/): Scripts for generating the genomic track plots for CDKN2B-AS1 (Supp. Figure 2).

Other relevant software related to the manuscript's analyses:
- [STARE](https://github.com/schulzlab/stare): Predict enhancer-gene interactions
- [SNEEP](https://github.com/SchulzLab/SNEEP/): Predict impact of SNVs on transcription factor binding sites
- [coloc](https://cran.r-project.org/web/packages/coloc/index.html): Test colocalisation of two genetic traits


