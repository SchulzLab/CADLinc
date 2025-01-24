import sys, os

if (len(sys.argv) < 4):
	print("python3 getSNPSNAP_input.py snps_hg19.txt, notMappedSNPs.txt,  outputfile ")
else:

	snpFile = sys.argv[1] # GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg19_16_02_lead_proxySNPs_swapped_uniq.txt
	notMappedFile = sys.argv[2]
	outputFile = sys.argv[3] # 


	notMapped = []
	with open(notMappedFile, 'r') as i:
		for line in i:
			line = line.strip().split('\t')
			notMapped.append(line[0][3:] + ":" + line[2])

	with open(snpFile, 'r') as i, open(outputFile, 'w') as o:
		o.write("lead_snp\n")
		for line in i:
			line = line.strip().split('\t')
			snp = line[0][3:] + ":" + line[2]
			if not snp in notMapped:
				o.write(snp + '\n')

