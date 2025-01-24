import sys, os


def cellTypeCheck(snps, chr_, subfolders):

	print(chr_)

	for elem in subfolders:
		counter = 0
		counter_all = 0
		celltype = elem.split("_")[-1]
		with open(elem + "/allRandomSNPsPos_hg19.txt", 'r') as i:
			if chr_ == "1":
				with open(elem + "/allRandomSNPs_withAlleleInfo_hg19.txt", 'w') as o:
					o.write("chr:end_hg19\tmajor\tminor\trsID\tMAF\tleadSNP_hg19\n")
					for line in i:
						line = line.strip().split('\t')
						if line[1].split(":")[0] == chr_: # matched snp of the current chr_
							counter_all += 1
							pos_hg19 = "chr" + line[1]
							## check if snp is in the 1000 genome project
							if pos_hg19 in snps:
								o.write(pos_hg19 + '\t' + snps[pos_hg19] + "\t" + line[0]  + '\n' )
							else:
								counter +=1
			else:
				with open(elem + "/allRandomSNPs_withAlleleInfo_hg19.txt", 'a') as o:
					for line in i:
						line = line.strip().split('\t')
						if line[1].split(":")[0] == chr_: # matched snp of the current chr_
							counter_all += 1
							pos_hg19 = "chr" + line[1]
							## check if snp is in the 1000 genome project
							if pos_hg19 in snps:
								o.write(pos_hg19 + '\t' + snps[pos_hg19] + "\t" + line[0]  + '\n')
							else:
								counter +=1
		print(celltype + " " + str(counter_all) + " " + str(counter))

if (len(sys.argv) < 3):
	print("python3 getRandomSNPs_SNEEP.py, snipaFile, randomSNPDir")
else:

	snipaFile = sys.argv[1] ## snipa file containing each snp once without LD structure /projects/sneep/work/SNEEP/GWAS_Catalog/snipa_data/snipa_snps_with_allele_info_sorted.txt
	randomSNPdir = sys.argv[2] ## result folder of the random snps from snpsnap and sneep
	# get the subfolders of the inputfolder
	subfolders = [ f.path for f in os.scandir(randomSNPdir) if f.is_dir() ]

	chr_ = "chr1"
	print(chr_)
	snps = {} 
	with open(snipaFile, 'r') as i: # read in snps from chr1
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			if line[0] != chr_:
				print("number snps: " + str(len(snps)))
				cellTypeCheck(snps, chr_[3:], subfolders) ## get allele freq for the cell type specific snps
				chr_ = line[0] ## adapte chr
				print(chr_)
				snps = {} # remove all elements from snps dictionary and add new entry
				snps[line[0] + ":" + line[2]] = line[10] + '\t' + line[8] + '\t' + line[6] + '\t' + line[9]  # hg19 pos chr:pos (1-based) plus the necessary info  chr:end major minor rsID MAF
			else:
				## add a snps to the list of the current chromosome
				snps[line[0] + ":" + line[2]] = line[10] + '\t' + line[8] + '\t' + line[6] + '\t' + line[9]  # hg19 pos chr:pos (1-based) plus the necessary info
			
	cellTypeCheck(snps, chr_[3:], subfolders) ## get allele freq for the cell type specific snps
