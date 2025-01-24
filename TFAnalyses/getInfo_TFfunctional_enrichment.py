import sys, os


TFs = ["NFIX","THRB", "CTCF", "PATZ1","TCFL5","KLF16", "TFEB", "ZNF610", "TFAP2E", "CTCFL"]
mapping = {"CTCF": ["CTCF_MA0139.1"], "THRB" : ["THRB_MA1575.1", "THRB_MA1574.1"], "NFIX" : ["NFIX_MA1528.1"]}

## add target gene information to heatmap data and combine the results for the motifs of TF THRB 


def readFile(filename, o, tf):
	with open(filename , 'r') as i:
		i.readline() # skip header
		for line in i:
			line = line.strip().split("\t")
			if line[0] == "GO:BP" or line[0] == "WP" or line[0] == "REAC" or line[0] == "GO:MF":
				o.write(tf + '\t' + line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[16] + '\n')


if len(sys.argv) < 3:
	print("python3 getTargetGenes.py pathToFiles  outputfile")
else:
	pathToData = sys.argv[1]
	outputFile = sys.argv[2]

	with open(outputFile , "w") as o:
		o.write("TF\tsource\tnative\tname\tp_value\tfraction\n")
		for elem in TFs:
			print(elem)

			if elem in mapping.keys():
				print("here")
				helper = mapping[elem]
				for h in helper:
					filename = pathToData + "FixDC_" + h + "_gProfilerDF.tsv" 
					readFile(filename, o, elem)

			else:
				filename = pathToData + "FixDC_" + elem + "_gProfilerDF.tsv" 
				readFile(filename, o, elem)


