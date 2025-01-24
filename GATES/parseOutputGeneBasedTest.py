import sys, os
import HelperFunctions

if (len(sys.argv) < 4):
	print("python3 parseOutputGeneBasedTest.py inputFile, outputFile")
else:
	inputFile = sys.argv[1] 
	annotationFile = sys.argv[2]
	outputFile = sys.argv[3]

	annotation = HelperFunctions.gene_biotypes(annotationFile)

	with open(inputFile, 'r') as i, open(outputFile, 'w') as o:
		i.readline() # skip first line
		o.write("GENE\tp.value\tNo_SNPs\tp.adj\tgene_type\tgene_name\n")
		for line in i:
			line = line.strip().split(',')

			if float(line[4]) <= 0.05:
				o.write(line[1][1:-1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' +  annotation[line[1][1:-1]]['gtf'] + '\t' + annotation[line[1][1:-1]]['name'] + '\n')

