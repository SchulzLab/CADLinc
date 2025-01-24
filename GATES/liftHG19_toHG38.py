import sys, os
from liftover import get_lifter

if (len(sys.argv) < 3):
	print("python3 liftHg19_hg38.py summaryStatisticNotFiltered.txt , outputFile")
else:
	inputFile = sys.argv[1]
	outputFile = sys.argv[2]

	converter = get_lifter('hg19', 'hg38')
	
	counter_lost = 0
	with open(inputFile, 'r') as i, open(outputFile, 'w') as o:
		header = i.readline() #skip header
		o.write(header)
		for line in i:
			line = line.strip().split('\t')
			chrom = line[1] 
			pos_start = int(line[2])-1
			pos_end = int(line[2])
			alleles = line[0].split("_") #1:61743	C	G
			con_start = converter[chrom][pos_start]
			con_end = converter[chrom][pos_end]
			if con_start == [] or con_end == [] or con_start[0][0] != con_end[0][0] or con_start[0][2] != con_end[0][2]:
				counter_lost += 1
			else:
				if con_start[0][2] == "+":
					gene_marker = con_end[0][0] + ":" + str(con_end[0][1]) + "_" + alleles[1] +   "_" + alleles[2]
					o.write(gene_marker + '\t' + con_end[0][0][3:] + '\t' + str(con_end[0][1])  + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] + '\t' + line[9] + '\t' + line[10] + '\t' + line[11] + '\t' + line[12] + '\t' + line[13] + '\t' + line[14] + '\t' + line[15] + '\t' + line[16] + '\t' + line[17] + '\t' + line[18] + '\t' + line[19] + '\t' + line[20] + '\n' )
				else:
					gene_marker = con_end[0][0] + ":" + str(con_end[0][1]+1) + "_" + alleles[1] +   "_" + alleles[2]
					o.write(gene_marker + '\t' + con_end[0][0][3:] + '\t' + str(con_end[0][1] + 1)  + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] + '\t' + line[9] + '\t' + line[10] + '\t' + line[11] + '\t' + line[12] + '\t' + line[13] + '\t' + line[14] + '\t' + line[15] + '\t' + line[16] + '\t' + line[17] + '\t' + line[18] + '\t' + line[19] + '\t' + line[20] + '\n' )

	print("lost snps: " + str(counter_lost))
