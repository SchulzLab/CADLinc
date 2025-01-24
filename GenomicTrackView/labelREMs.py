import sys, os


#GENE="ENSG00000259673"


def labelFile(file_, c,o):

	with open(file_, 'r') as i:
		for line in i:
			flag = "-"
			line = line.strip().split('\t')	
			genes = line[4].split(",")
			if GENE in genes: 
				flag = "+"

			o.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + c + '\t' + flag + '\n' )
		

if (len(sys.argv) < 5):
	print("python3 labelREMs.py file1, file2, file3,  celltypes,  outputFile" ) 
else:
	files = sys.argv[1] ## comma seperated list of files
	celltypes = sys.argv[2] # comma seperated list of celltypes (no spaces), when considering all, define "ALL"
	GENE = sys.argv[3]
	outputFile = sys.argv[4]

	files = files.split(',')
	print(files)
	celltypes = celltypes.split(",")
	print(celltypes)	

	counter = 0
	with open(outputFile, 'w') as o:
		o.write("chr\tstart\tend\tid\tgene\tgroup\tstrand\n")
		for f in files: 
			labelFile(f, celltypes[counter], o)
			counter += 1
