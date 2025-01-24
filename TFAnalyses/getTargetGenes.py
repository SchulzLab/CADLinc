import sys
import HelperFunctions


def countGenes(genes, mappingGene):

	count = 0
	for t in genes:
		biot = mappingGene[t]['general']
		if biot == "non-coding RNA" or biot == "protein-coding":
			count +=1
	return count


## add target gene information to heatmap data and combine the results for the motifs of TF THRB 

if (len(sys.argv) < 4):
	print("python3 getTargetGenes.py inputFile, TF_table_extended, annotationFile, outputfile")
else:
	inputFile = sys.argv[1]
	TfTableFile = sys.argv[2]
	annotationFile = sys.argv[3]
	outputFile = sys.argv[4]


	## read in TF table
	known_CAD = {} # TF to number of known CAD loci
	new_CAD = {} # TF to number of new CAD loci (not including known CAD loci)
	targetGenes = {} # TF to number of target genes minus new_CAD and known_CAD
	numberRegSNVs  = {} # TF to number of regSNVs in celltypes with high odds ratio
	ncGenes = {} # TF to number of ncGenes
	cGenes = {} # TF to number of cGenes
	otherGenes = {} # TF to number of pseudo or TEC genes
	mappingGene = HelperFunctions.gene_biotypes(annotationFile) ## per gene the info whether it is coding or not
	helper_cad = set()
	helper_new_cad = set()
	helper_targetGenes = set()
	helper_regSNVs = set()
	helper_ncGenes = set()
	helper_cGenes = set()
	helper_others = set()

	targetGeneNames = ""

	overallKnownCAD_genes = set()
	overallCandidateCAD_genes = set()

	interestingTFs = []
	with open(inputFile, 'r') as i:
		header = i.readline() 
		header = header.split(",")[2:]
		print(header) 
		interestingTFs = header

	with open(TfTableFile, 'r') as i:
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			tf = line[0].split("(")[0] ## remove (ID info)

			## initialize ncGenes and cGenes for the current TF
			if not tf in ncGenes.keys():
				ncGenes[tf] = 0
			if not tf in cGenes.keys():
				cGenes[tf] = 0
			if not tf in otherGenes.keys():
				otherGenes[tf] = 0

			if tf == "THRB": ## combine info of this TF (two TF motifs with high odds ratio)
				helper_targetGenes.update(set(line[11].split(",")))
				helper_cad.update(set(line[15].split(",")[:-1]))
				helper_new_cad.update(set(line[16].split(",")[:-1]))
				helper_regSNVs.update(set(line[13].split(",")))
	
				targetGeneNames = line[11].split(",")
				for t in targetGeneNames: 
					biot = mappingGene[t]['general']
					if biot == "non-coding RNA":
						helper_ncGenes.add(t)
					if biot == "protein-coding":
						helper_cGenes.add(t)
				
			else:
				numberRegSNVs[tf] = len(line[13].split(","))
				known_CAD[tf] = countGenes(line[15].split(",")[:-1], mappingGene)
				new_CAD[tf] = countGenes(line[16].split(",")[:-1], mappingGene)- known_CAD[tf]
				targetGenes[tf] = countGenes(line[11].split(","), mappingGene) - known_CAD[tf] - new_CAD[tf]
				if tf in interestingTFs: 
					overallKnownCAD_genes = overallKnownCAD_genes.union(set(line[15].split(",")[:-1]))
					overallCandidateCAD_genes = overallCandidateCAD_genes.union(set(line[16].split(",")[:-1]))
					
				targetGeneNames = line[11].split(",")
				for t in targetGeneNames: 
					biot = mappingGene[t]['general']
					if biot == "non-coding RNA":
						ncGenes[tf] += 1
					if biot == "protein-coding":
						cGenes[tf] += 1
					
					#else:
					#	otherGenes[tf] += 1
						#print("TF: " + tf + " biotype: " + biot + " gene: " + t)

	## add info for THRB
	numberRegSNVs["THRB"] = len(helper_regSNVs)
	known_CAD["THRB"] = countGenes(helper_cad, mappingGene)
	new_CAD["THRB"] = countGenes(helper_new_cad - helper_cad, mappingGene)
	targetGenes["THRB"] = countGenes(helper_targetGenes, mappingGene) - new_CAD["THRB"] - known_CAD["THRB"] 
	overallKnownCAD_genes = overallKnownCAD_genes.union(set(helper_cad))
	overallCandidateCAD_genes = overallCandidateCAD_genes.union(set(helper_new_cad))

	print(str(len(overallKnownCAD_genes)) + " " + str(overallKnownCAD_genes))
	print(str(len(overallCandidateCAD_genes)) + " " + str(overallCandidateCAD_genes))

	cGenes["THRB"] = len(helper_cGenes)
	ncGenes["THRB"] = len(helper_ncGenes)
	otherGenes["THRB"] = len(helper_others)
			
	with open(inputFile, 'r' ) as i, open(outputFile, 'w') as o:

		header = i.readline().strip().split(",")
		## write new header and exclude THRB at position 3 and 22 
		for elem in header[:3]:
			o.write(elem + ",") 
		for elem in header[4:22]:
			o.write(elem + ",") 
		for elem in header[23:]:
			o.write(elem + ",") 
		o.write("THRB\n") ## add new entry for combined info

		tf_order = header[2:] 

		## add additional TF specific information as rows at the beginning of the new file
		row_known_CAD = "knownCAD,knownCAD" ## celltype info and label
		row_new_CAD = "newCAD,newCAD" ## celltype info and label
		row_targetGenes = "targetGenes,targetGenes" ## celltype info and label
		row_regSNVs = "regSNVs,regSNVs" ## celltype info and label
		row_ncGenes = "ncGenes,ncGenes"
		row_cGenes = "cGenes,cGenes"
		row_others = "others,others"
		for elem in tf_order: 
			#print(elem)
			if elem == "THRB(MA1575.1)" or elem == "THRB(MA1574.1)": ## combine THRB
				continue
			else:
				row_known_CAD += "," + str(known_CAD[elem])
				row_new_CAD += "," + str(new_CAD[elem])
				row_targetGenes += "," + str(targetGenes[elem])
				row_regSNVs += "," + str(numberRegSNVs[elem])
				row_ncGenes += "," + str(ncGenes[elem])
				row_cGenes += "," + str(cGenes[elem])
				row_others += "," + str(otherGenes[elem])

				##check if counts are correct
				if (known_CAD[elem] + new_CAD[elem] + targetGenes[elem]) != (ncGenes[elem] + cGenes[elem] + otherGenes[elem]):
					print("AHHHH: " + elem)
		
		## add THRB
		row_known_CAD += "," + str(known_CAD["THRB"])
		row_new_CAD += "," + str(new_CAD["THRB"])
		row_targetGenes += "," + str(targetGenes["THRB"])
		row_regSNVs += "," + str(numberRegSNVs["THRB"])
		row_ncGenes += "," + str(ncGenes["THRB"])
		row_others += "," + str(otherGenes["THRB"])
		row_cGenes += "," + str(cGenes["THRB"])
		##check if counts are correct
		if (known_CAD["THRB"] + new_CAD["THRB"] + targetGenes["THRB"]) != (ncGenes["THRB"] + cGenes["THRB"] + otherGenes["THRB"]):
			print("AHHHH: THRB")
		
		o.write(row_known_CAD + '\n' + row_new_CAD  + '\n' + row_targetGenes + '\n' + row_regSNVs + '\n' + row_ncGenes + '\n' + row_cGenes + '\n' + row_others + '\n')
		
		## keep the info from the inputFile in the outputFile
		for line in i:

			## write entries but not THRB
			line = line.strip().split(",")
			for elem in line[:3]:
				o.write(elem + ",") 
			for elem in line[4:22]:
				o.write(elem + ",") 
			for elem in line[23:]:
				o.write(elem + ",") 

			#combine THRB info
			o.write(str(max(float(line[3]), float(line[22]))) + '\n')




