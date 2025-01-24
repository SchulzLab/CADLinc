import sys, os
import pickle


if (len(sys.argv) < 6):
	print("python3 check_eQTL_data.py TF_info, eQTLs_GTEx, eQTLs_STARNET,  outputFile" ) 
else:
	cellTypeSpecificSneepResult = sys.argv[1]  
	eQTLs_GTEx_file = sys.argv[2]
	eQTLs_STARNET_file = sys.argv[3]
	openGenes = sys.argv[4] # OpenGenesCell.pkl  
	outputFile = sys.argv[5]

	## read in  active genes per celltype
	activeGenes = pickle.load(open(openGenes, 'rb'))  ## celltype to set of open genes

	##read in celltype specific snp info to gene 

	mapping = {} # per celltype provide for each snp position a list with the associated genes
	subfolders = [ f.path for f in os.scandir(cellTypeSpecificSneepResult) if f.is_dir() ]
	for f in subfolders:
		c = f.split("_")[4]
		pos_genes = {} # for each position provide the asssociated genes and their gene names
		with open(f + "/result.txt", 'r') as i:
			
			for line in i:
				line = line.strip().split('\t')
				if not line[0] in pos_genes.keys():
					genes = line[15].split(",")
					gene_names = line[16].split(",")
					pos_genes[line[0]] = {'genes' : genes, 'gene_names' : gene_names}
		mapping[c] = pos_genes	

	## read in eQTLs from GTEx
	with open(eQTLs_GTEx_file, 'r') as i, open(outputFile, 'w') as o:	
		o.write("celltype\tTF\trsID\tpos\tgene\tgeneName\tsource\n")
		i.readline() # skip header
		for line in i:
			line = line.strip().split(",")
			c = line[1][1:-1]
			tf = line[2][1:-1]
			rsid = line[3][1:-1]
			pos = line[4][1:-1].split(":")
			pos = pos[0] + ":" + pos[1]
			gene = line[6][1:-1].split(".")[0]
			gene_name = line[7][1:-1]

			#check if the associated gene is the same as foru our prediction
			pos_genes = mapping[c]
			if pos in pos_genes.keys():
				if gene in pos_genes[pos]['genes']: #or gene_name in pos_genes[pos]['gene_names']:
					#check if the gene is active in the celltype
					if gene in activeGenes[c]:
						print("match" + str(line))

	## read in STARNET data 
	with open(eQTLs_STARNET_file, 'r') as i, open(outputFile, 'a') as o:
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			pos_genes = mapping[line[0]]

			pos = line[3].split(":")
			pos = pos[0] + ":" + pos[1]

			gene = line[9]

			if pos in pos_genes.keys():
				if gene in pos_genes[pos]['genes'] and gene in activeGenes[c]:
					o.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + "\t" + gene + '\t' + "-" + "\tSTARNET\n" )

