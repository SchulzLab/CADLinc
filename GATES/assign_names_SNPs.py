import sys, os
if (len(sys.argv) < 4):
	print("python3 assign_names_SNPs.py PLINK file, SNP file, outputFile ")
else:
	PLINKFile = sys.argv[1] 
	SNVFile = sys.argv[2]
	outputFile = sys.argv[3]

	plinkSNVs = set() # chr to pos
	plinkSNVs_info = {} #chr_pos to full info
	with open(PLINKFile, 'r' ) as i:
		i.readline() # skip header
		for line in i:
			l = line.strip().split('\t')
			plinkSNVs.add(l[0] + ":" + l[3])	
			plinkSNVs_info[l[0] + ":" + l[3]] = line.strip()

	counter =0
	with open(SNVFile, 'r') as i, open(outputFile, 'w') as o:
		o.write("CHR\tSNP_ID\t0\tBP\tREF\tALT\tMarkerName\tCHR\tBP\tAllele1\tAllele2\tFreq1\tFreqSE\tMinFreq\tMaxFreq\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tCases\tEffective_Cases\tN\tMeta_analysis\tGENE\n")
		i.readline() # skip header
		for line in i:
			l = line.strip().split('\t')
			if l[1] + ":" + l[2] in plinkSNVs:
				counter+=1 
				o.write(plinkSNVs_info[l[1] + ":" + l[2]]  + '\t' + line)

	print("counter: " + str(counter))


