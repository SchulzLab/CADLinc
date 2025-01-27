#! /bin/bash

regions=$1 
key=$2

for sample in "AD" "BFUE" "BM" "Bcell" "CD4Tcell" "CD8Tcell" "CFUE" "CLP" "CM" "CMP" "EBaso" "GMP" "HC" "HSC" "KEC" "LBaso" "LMPP" "MEP" "MPP" "MSC" "Mac" "Mono" "NK" "Neutro" "Ortho" "Poly" "ProE" "SKMC" "VSMC" "aCM" "adFB" "arEC" "arFB" "cAD" "cEC" "cFB" "cLC" "cMac" "cNR" "cSMC" "mDC" "pDC" "sccFB" "vCM"  ;do

	bedtools intersect -b ${regions} -a ../enhancer_interactions_avghg38_unique/${sample}_*.bed -wo > ../figure_tables/${key}_${sample}.txt

	echo " ../figure_tables/${key}_${sample}.txt"

done

python labelREMs.py ../figure_tables/${key}_AD.txt,../figure_tables/${key}_BFUE.txt,../figure_tables/${key}_BM.txt,../figure_tables/${key}_Bcell.txt,../figure_tables/${key}_CD4Tcell.txt,../figure_tables/${key}_CD8Tcell.txt,../figure_tables/${key}_CFUE.txt,../figure_tables/${key}_CLP.txt,../figure_tables/${key}_CM.txt,../figure_tables/${key}_CMP.txt,../figure_tables/${key}_EBaso.txt,../figure_tables/${key}_GMP.txt,../figure_tables/${key}_HC.txt,../figure_tables/${key}_HSC.txt,../figure_tables/${key}_KEC.txt,../figure_tables/${key}_LBaso.txt,../figure_tables/${key}_LMPP.txt,../figure_tables/${key}_MEP.txt,../figure_tables/${key}_MPP.txt,../figure_tables/${key}_MSC.txt,../figure_tables/${key}_Mac.txt,../figure_tables/${key}_Mono.txt,../figure_tables/${key}_NK.txt,../figure_tables/${key}_Neutro.txt,../figure_tables/${key}_Ortho.txt,../figure_tables/${key}_Poly.txt,../figure_tables/${key}_ProE.txt,../figure_tables/${key}_SKMC.txt,../figure_tables/${key}_VSMC.txt,../figure_tables/${key}_aCM.txt,../figure_tables/${key}_adFB.txt,../figure_tables/${key}_arEC.txt,../figure_tables/${key}_arFB.txt,../figure_tables/${key}_cAD.txt,../figure_tables/${key}_cEC.txt,../figure_tables/${key}_cFB.txt,../figure_tables/${key}_cLC.txt,../figure_tables/${key}_cMac.txt,../figure_tables/${key}_cNR.txt,../figure_tables/${key}_cSMC.txt,../figure_tables/${key}_mDC.txt,../figure_tables/${key}_pDC.txt,../figure_tables/${key}_sccFB.txt,../figure_tables/${key}_vCM.txt  AD,BFUE,BM,Bcell,CD4Tcell,CD8Tcell,CFUE,CLP,CM,CMP,EBaso,GMP,HC,HSC,KEC,LBaso,LMPP,MEP,MPP,MSC,Mac,Mono,NK,Neutro,Ortho,Poly,ProE,SKMC,VSMC,aCM,adFB,arEC,arFB,cAD,cEC,cFB,cLC,cMac,cNR,cSMC,mDC,pDC,sccFB,vCM ENSG00000240498  ../figure_tables/rSNVs_CADlinc_${key}_allCelltypes.txt

