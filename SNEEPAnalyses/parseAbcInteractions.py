import sys, os


if (len(sys.argv) < 3):
	print("python3 parseAbcInteractions.py pathToInteractionsFile, outputFile ")
else:

	interactionDir = sys.argv[1]
	outputFile = sys.argv[2]

	with open(interactionDir + "/celltypes.txt" , 'r') as c, open(outputFile, 'w') as o:

		#iterate over all files
		for line in c: 
			line = line.strip()
			celltype = line.strip().split("_")[0]
			print(celltype)

			# open file per celltype
			with open(interactionDir + "/" + line, 'r') as i:
				
				for line2 in i:

					line2 = line2.strip().split('\t')
					if line2[4] != "-":
						genes = line2[4].split(",")
						for elem in genes:
							o.write(line2[0] + "\t" + line2[1] + "\t" + line2[2] + "\t" + elem + "\t" + celltype + "\t-\t-\t-\t-\t-\t-\t-\n")
