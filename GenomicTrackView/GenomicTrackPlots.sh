#! /bin/bash

##parsing script that produces an outputFile only including the REMs for IQCH-AS1 in the format chr	start	end	id	gene	group, where group is the used celltype

## removed the LD SNVs from the file ../figure_tables/rSNVs_CDKN2B_AS1.txt -> figure_tables/lead_TF_SNVs_CDKN2B_AS1.txt
bash workflow_intersect.sh ../figure_tables/lead_TF_SNVs_CDKN2B_AS1.txt CDKN2B_AS1_leadSNVs_all_1312 # output file rSNVs_CADlinc_CDKN2B_AS1_leadSNVs_all_1312_allCelltypes.txt
python labelREMs.py ../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_CD8Tcell.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_CD4Tcell.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_EBaso.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_LBaso.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_mDC.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_ProE.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_GMP.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_CMP.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_MPP.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_vCM.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_aCM.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_cAD.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_cNR.txt,../figure_tables/CDKN2B_AS1_leadSNVs_all_1312_CM.txt  CD8Tcell,CD4Tcell,EBaso,LBaso,mDC,ProE,GMP,CMP,MPP,vCM,aCM,cAD,nNR,CM   ENSG00000240498  ../figure_tables/rSNVs_CADlinc_CDKN2B_AS1_leadSNVs_all_1312_allCelltypes.txt

Rscript genomicView_CDKN2B_AS1_2.0.R  ../figure_tables/rSNVs_CADlinc_CDKN2B_AS1_leadSNVs_all_1312_allCelltypes.txt  ../figure_tables/lead_TF_SNVs_CDKN2B_AS1.txt ../base_data/gencode.v38.annotIation.gtf genomicView_CDKN2B-AS1_leadSNVs_13_12.pdf

