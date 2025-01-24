# create txt files with groups of SNPs based on their genes
## save them in Files_Gene_Based

import sys, os

if (len(sys.argv) < 4):
	print("python3 group_SNP_genes_Enhancers.py SNP file, outputFile, OUtputgeneFile ")
else:
	SNVFile = sys.argv[1]
	outputDir = sys.argv[2]
	outputGeneFile = sys.argv[3]

	g = ""
	helper = set() #store rsID per gene
	overallGenes = set()
	counterGenes = 0
	counterREMs = 0
	counter = 0
	chrom = 0
	with open(SNVFile, 'r') as i:
		for line in i:
			line = line.strip().split('\t')
			if line[0] != "CHR": # is not header
				gene = line[-1]
				if g == gene or g == "": 
					if not line[1] in helper:
						helper.add(line[1]) ## ad rsID
					else: 
						counter +=1
					g = gene
					chrom = line[0] #set chrom
				else:
					counterGenes +=1
					overallGenes.add(g)
					## write info from old gene
					with open(outputDir + "/" + chrom + "/" + g + ".txt", 'w') as o:
						for elem in helper:
							counterREMs += 1
							o.write(elem + '\n')
					## reset for new gene
					helper = set()
					helper.add(line[1])
					g = gene
					chrom = line[0]
				
	with open(outputDir + "/" + chrom + "/" + g + ".txt", 'w') as o:
		for elem in helper:
			counterREMs += 1
			o.write(elem + '\n')
	counterGenes +=1
	overallGenes.add(g)

	with open(outputGeneFile, 'w') as o:
		o.write("gene\n")
		for elem in overallGenes:
			o.write(elem + '\n') 
	print("number genes: " + str(counterGenes) +  " and length list of genes: " + str(len(overallGenes)))
	print("counter SNPs considered: " + str(counterREMs))	
	print("counter SNPs duplicated rsID but differ in alleles: " + str(counter))

