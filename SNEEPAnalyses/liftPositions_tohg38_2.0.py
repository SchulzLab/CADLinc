import sys, os
from liftover import get_lifter

### SNP positions in input file are 0-baseds
## pyllift need 1 - based posiyions

##lift uniq lead an proxy SNPs to hg38

if (len(sys.argv) < 3):
	print("python3 liftPositions_tohg38.py leadAndProxySNPs.txt , outputFile")
else:
	leadSNPsFile = sys.argv[1]
	outputFile = sys.argv[2]

	converter = get_lifter('hg19', 'hg38')

	#read in file
	with open(leadSNPsFile, 'r') as i, open(outputFile, 'w') as o:
		o.write("#chr\tstart\tend\tallele1\tallele2\trsId\tMAF\trsIDsLead\tstatus\tposition_hg19\n")

		for line in i:
			line = line.strip().split('\t') # chr start end a1 a2 rsID MAF rsIDLead status 
			chrom = line[0][3:] #remove 'chr'
			position_end = int(line[2])
			position_start = int(line[1])
			hg19_pos = line[0] + ":" + line[1] + "-" + line[2]
			con_start = converter[chrom][position_start]
			con_end = converter[chrom][position_end]
			## check if the converter worked, if the positions are on the same chromosome and if the strand info is the same
			if con_start == [] or con_end == [] or con_start[0][0] != con_end[0][0] or con_start[0][2] != con_end[0][2]:
				print(line)
			else:
				if con_start[0][2] == "+":
					o.write(con_start[0][0] + '\t' + str(con_start[0][1]) + '\t' + str(con_end[0][1]) + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] +  '\t' + hg19_pos + '\n')
				else:
					o.write(con_start[0][0] + '\t' + str(con_end[0][1] + 1) + '\t' + str(con_start[0][1] + 1) + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + line[8] +  '\t' + hg19_pos + '\n')
					
