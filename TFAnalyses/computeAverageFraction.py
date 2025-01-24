import sys, os

if (len(sys.argv) < 3):
	print("python3 getTargetGenes.py inputFile  outputfile")
else:
	inputFile = sys.argv[1] # TFs_functionalEnrichment_detailedPathways.txt 
	outputFile = sys.argv[2]

	TF = "NFIX"
	pathway = "Immune response and inflammation"
	source = "GO:BP"
	counter = 0
	average_p = 0
	average_f = 0
	with open(inputFile , "r") as i, open(outputFile, 'w') as o:
		o.write("TF\tsource\tenrichedPathway\taveragePvalue\taverageFraction\n")
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			current_TF = line[0]
			current_pathway = line[6]
			if current_TF == TF:
				if current_pathway == pathway: 
					average_p += float(line[4])
					average_f += float(line[5])
					counter += 1
					source = line[1]
				else: ## compute average when pathway differes inbetween a TF 
					o.write(TF + '\t' + source + '\t' + pathway + '\t' + str(average_p/counter) + '\t' + str(average_f/counter) + '\n')
					counter = 1
					source = line[1]
					pathway = current_pathway
					average_p = float(line[4])
					average_f = float(line[5])

			else: ## compute average
				o.write(TF + '\t' + source + '\t' + pathway + '\t' + str(average_p/counter) + '\t' + str(average_f/counter) + '\n')
				counter = 1
				source = line[1]
				TF = current_TF
				pathway = current_pathway
				average_p = float(line[4])
				average_f = float(line[5])

		## write last entry
		o.write(TF + '\t' + source + '\t' + pathway + '\t' + str(average_p/counter) + '\t' + str(average_f/counter) + '\n')


