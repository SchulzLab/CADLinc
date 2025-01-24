import sys, os

if (len(sys.argv) < 3):
	print("python3 parseSNIPA_SNPs.py snipaDir, outputfile")
else:

	snipaDir = sys.argv[1]
	outputFile = sys.argv[2]

	with open(outputFile, 'w') as o:
		for i in range(1,24):
			print(i)
			filename = ""
			if (i == 23):
				filename = "grch37-1kgpp3v5-eur-chrX-ld"
			else:
				filename = "grch37-1kgpp3v5-eur-chr" + str(i) + "-ld"

			with open(snipaDir + "/" + filename, 'r') as i:
				if i == 1:	
					o.write(i.readline())
				else:
					i.readline() # skip header
				for line in i:
					o_line = line
					line =  line.strip().split('\t')
					if line[1] == line[2]:
						o.write(o_line)


