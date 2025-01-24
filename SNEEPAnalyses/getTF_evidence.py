import sys, os

def checkOverlap(tf_genes, geneList):

	overlap_tf_genes = {} # per TF the overlapping genes
	for tf in tf_genes.keys():
		collect = set()
		g = tf_genes[tf]
		for elem in g: 
			if elem in geneList:
				collect.add(elem)
		overlap_tf_genes[tf] = collect
	return overlap_tf_genes

def checkElem(o, geneList, elem):

	if elem in geneList:
		o.write('\t' + "True")
	else:
		o.write('\t' + "False")
	
if (len(sys.argv) < 9):
	print("python3 getTF_evidence.py TF_table_newestVersion, TFevidenceDir, TF mapping") 
else:
	tfTableFile = sys.argv[1]  
	tfEvidenceDir = sys.argv[2]
	mappingFileTF_ensemble = sys.argv[3]
	coExpressedGenesFile = sys.argv[4]
	ourCADLociFile = sys.argv[5]
	outputFileCAD = sys.argv[6]
	outputFileTFs = sys.argv[7]
	outputFileTFTable = sys.argv[8]

	## read in TF mapping	
	TF_mapping = {} # TF name to id
	with open(mappingFileTF_ensemble, 'r') as i:
		i.readline()
		for line in i:
			line = line.strip().split('\t')
			TF_mapping[line[1]] = line[0]

	## read in information from TF table
	tfs = set() # TF names removed () and splitted ::
	tfs_id = set() # TF ensembl ids
	tf_genes = {} # per TF get all target genes 

	tf_table = {} # TF to info from table
	#genes_mapping = {} # ensembl ID ti name mapping
	with open(tfTableFile, 'r') as i: 
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			o_tf = line[0]
			tf_table[o_tf] = line[1:]
			tf = line[0].split("(")[0]
			if "::" in tf:
				tf = tf.split("::")
				for elem in tf:
					tfs.add(elem)
					tfs_id.add(TF_mapping[elem])

			else:
				tfs.add(tf)
				tfs_id.add(TF_mapping[tf])
			# get gene info 
			tf_genes[o_tf] = line[7].split(",")

	## read in coexpression data: 
	coExpressedGenes_known = []
	coExpressedGenes_new = []
	with open(coExpressedGenesFile, 'r') as i:
		i.readline() #skip header
		for line in i:
			line = line.strip().split('\t')
			if line[5] != "-" and float(line[5]) <= 0.05:
				if line[2] == "False":
				
					coExpressedGenes_new.append(line[0])
				else:
					coExpressedGenes_known.append(line[0])
		
	overlap_coexpressedGenes_known = checkOverlap(tf_genes, coExpressedGenes_known)
	overlap_coexpressedGenes_new = checkOverlap(tf_genes, coExpressedGenes_new)

	## read in our predicted CAD loci 
	ourCADGenes = set()
	with open(ourCADLociFile, 'r') as i: 
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			ourCADGenes.add(line[0])
	
	overlap_ourCADGenes = checkOverlap(tf_genes, list(ourCADGenes))

	## read in CAD genes from zhifen
	files = os.listdir(tfEvidenceDir + "/CAD_genes_Zhifen//")
	#print(files)
	GWAS = set()
	TWAS = set()
	TWAS_GWAS = set()
	for f in files:	
		with open(tfEvidenceDir + "/CAD_genes_Zhifen/" + f, 'r') as i:
			print(f)
			for line in i:
				line = line.strip()
				if "GWAS" in f and not "TWAS" in f:
					GWAS.add(line)
				elif "TWAS" in f and not "GWAS" in f:
					TWAS.add(line)
				else:
					TWAS_GWAS.add(line)

	print(len(GWAS))
	print(len(TWAS))
	print(len(TWAS_GWAS))
	overlap_GWAS = checkOverlap(tf_genes, GWAS)
	overlap_TWAS = checkOverlap(tf_genes, TWAS)
	overlap_TWAS_GWAS = checkOverlap(tf_genes, TWAS_GWAS)

	overlap_CAD = {} #TF to overlapped genes 	
	with open(outputFileCAD, 'w') as o:
		o.write("TF\tgene\tGWAS\tTWAS\tGWAS_TWAS\n")
		for tf in tf_genes.keys():
			seen = set()
			gwas = overlap_GWAS[tf]
			twas = overlap_TWAS[tf]
			both = overlap_TWAS_GWAS[tf]

			for elem in both:
				seen.add(elem)
				o.write(tf + '\t' + elem) 
				checkElem(o, gwas, elem)
				checkElem(o, twas, elem)
				o.write("\tTrue\n")

			overlap_CAD[tf] = list(seen)

			for elem in gwas:
				if not elem in seen:
					seen.add(elem)
					o.write(tf + '\t' + elem + "\tTrue")
					checkElem(o, twas, elem)
					checkElem(o, both, elem)
					o.write('\n')

			for elem in twas:
				if not elem in seen:
					seen.add(elem)
					o.write(tf + '\t' + elem + "\tFalse\tTrue\tFalse\n")


	## comparison to single cell immune landscape data 
	## check if one of our TFs is differential expressed in blood vs plaque immune cells in atherosclerosis

	files = os.listdir(tfEvidenceDir + "/SupplementaryData_scImmuneLandscape_41591_2019_590_MOESM3_ESM/")
	print(files)
	blood_vs_plaque = set()
	patients_ASYM_vs_SYM = set()
	for f in files:
		helper = f.split("_")[0]
		with open(tfEvidenceDir + "/SupplementaryData_scImmuneLandscape_41591_2019_590_MOESM3_ESM/" + f, 'r') as i:
			## skip the first four lines
			i.readline() 
			i.readline()
			i.readline()
			i.readline()
			for line in i:
				line = line.strip().split(";")
				if line[0] in tfs:
					if "6" in helper:
						patients_ASYM_vs_SYM.add(line[0])
					else:
						blood_vs_plaque.add(line[0])

	print("blood vs plaque: " + str(blood_vs_plaque))
	print("patient data: " + str(patients_ASYM_vs_SYM))

	with open(outputFileTFs, 'w') as o:
		o.write("TF\tbloodVsPlaque\tASYMvsSYM\n")
		seen = set()
		for elem in blood_vs_plaque:
			o.write(elem + '\t' + "True")
			checkElem(o, patients_ASYM_vs_SYM, elem)
			o.write("\n")
		for elem in patients_ASYM_vs_SYM:
			if not elem in seen:
				seen.add(elem)
				o.write(elem + '\t' + "\tFalse\tTrue\n")
				
	with open(outputFileTFTable, 'w') as o:
		o.write("TF\topenPromoter\tTF_singleCellImmuneLandscape\t#regSNVsOverallCelltypes\tregSNVsOverCellTypesHighOddsRatio\tcelltypesOddsRatioHigher5MoreThan5RegSNVs\t#targetGenes\t#targetGenesCADAssociated\t#targetGenesNewCADGenes\t#targetGenesCo_expr_knownLoci\t#targetGenesCo_expr_newLoci\ttargetGenes\ttargetGenesName\trsIDs\tSNVsDetails\ttaregtGenesCadAssociated\ttargetGenesNewCadGenes\ttargetGenesCo_expr_knownLoci\t#targetGenesCo_expr_newLoci\n")
		for tf in tf_table:
			o.write(tf + "\t")
			info_ =  tf_table[tf]
			o.write(info_[0] + '\t')

			## write TFs single cell 
			tf_name = tf.split("(")[0]
			if "::" in tf_name:
				tf_name = tf_name.strip().split("::")
				print(tf_name)
				if tf_name[0] in blood_vs_plaque and tf_name[1] in blood_vs_plaque:
					
					if tf_name[0] in  patients_ASYM_vs_SYM and tf_name[1] in  patients_ASYM_vs_SYM:
						o.write("both\t")
					else:
						o.write("bloodVsPlaque\t")
				elif tf_name[0] in  patients_ASYM_vs_SYM and  tf_name[1] in  patients_ASYM_vs_SYM:
					o.write("ASYMvsSYM\t")
				else:
					o.write("-\t")
				
			else:
				if tf_name in blood_vs_plaque and tf_name in  patients_ASYM_vs_SYM:
					o.write("both\t")
				elif tf_name in blood_vs_plaque and not tf_name in  patients_ASYM_vs_SYM:
					o.write("bloodVsPlaque\t")
				elif not tf_name in blood_vs_plaque and tf_name in  patients_ASYM_vs_SYM:
					o.write("ASYMvsSYM\t")
				else:
					o.write("-\t")

			for elem in info_[1:4] :
				o.write(elem + '\t')
			
			o.write(info_[5] + '\t')

			o.write(str(len(overlap_CAD[tf])) + '\t'+ str(len( overlap_ourCADGenes[tf]) - len(overlap_CAD[tf])) + '\t' + str(len(overlap_coexpressedGenes_known[tf])) + '\t' +str(len(overlap_coexpressedGenes_new[tf])) + '\t' )

			for elem in info_[6:]:
				o.write(elem + '\t')

			## write number of CAD genes and which
			o_genes = ""
			for elem in overlap_CAD[tf]:
				o_genes += elem + ","
	

			## write target genes of our cad genes
			c_genes = "" 
			for elem in  overlap_ourCADGenes[tf]:
				c_genes += elem + ","
			

			#write co-expressed genes  new
			co_genes_new = ""
			for elem in  overlap_coexpressedGenes_new[tf]:
				co_genes_new += elem + ","

			#write co-expressed genes known
			co_genes_known = ""
			for elem in  overlap_coexpressedGenes_known[tf]:
				if not elem in overlap_ourCADGenes[tf]:
					co_genes_known += elem + ","

			o.write(o_genes +'\t' + c_genes + '\t' + co_genes_known + '\t' + co_genes_new +   '\n')


			
	
