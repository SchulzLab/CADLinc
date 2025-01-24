import sys
import HelperFunctions


if len(sys.argv) < 4:
	print("python3 getNumberncGenes_Fig3E.py inputFile annotationFile,  outputfile")
else:
	inputFile = sys.argv[1]
	annotationFile = sys.argv[2]
	outputFile = sys.argv[3]
	
	mappingGene = HelperFunctions.gene_biotypes(annotationFile)  # per gene the info whether it is coding or not

	genes = set()
	counter_nc = 0
	counter_c = 0
	counter_others = 0
	with open(inputFile, 'r') as i: 
		i.readline()  # skip header
		for line in i:
			line = line.strip().split('\t')
			g = line[2]
			if not g in genes: 
				genes.add(g)
				
				biot = mappingGene[g]['general']
				if biot == "non-coding RNA":
					counter_nc += 1
				elif biot == "protein-coding":
					counter_c += 1 
				else:
					counter_others += 1

	with open(outputFile, 'w') as o:
		o.write("Biotype\tcount\nnon-coding RNA\t" + str(counter_nc) + '\n' + "protein-coding\t" + str(counter_c) + "\nothers\t" + str(counter_others) + '\n')
	




	







