import sys, os
from liftover import get_lifter
import random 
from random import randrange

BACKGROUND_ROUNDS = 500


def lift_hg19_to_hg38(converter, inputFile, outputFile):

	hg38_pos = {}
	counter_lost = 0
	with open(inputFile, 'r') as i, open(outputFile, 'w') as o:
		o.write("chr\tstart\tend\tmajor\tminor\trsID\tMAF\tleadSNP_hg19\n")
		i.readline() #skip header
		for line in i:
			line = line.strip().split('\t')
			chrom = line[0].split(":")[0][3:] #remove 'chr'
			pos_end = int(line[0].split(":")[1]) ## 1-based position
			pos_start = int(line[0].split(":")[1])-1 ## 0-based position
			con_start = converter[chrom][pos_start]
			con_end = converter[chrom][pos_end]
			if con_start == [] or con_end == [] or con_start[0][0] != con_end[0][0] or con_start[0][2] != con_end[0][2]:
				counter_lost += 1
			else:
				if con_start[0][2] == "+":
					new_pos = con_start[0][0] + '\t' + str(con_start[0][1]) + '\t' + str(con_end[0][1])
					new_pos_l = [con_start[0][0], str(con_start[0][1]), str(con_end[0][1])]
				else: ## for the negative strand the positions need to be shifted by plus one
					new_pos = con_start[0][0] + '\t' + str(con_end[0][1] + 1) + '\t' + str(con_start[0][1] + 1)
					new_pos_l = [con_start[0][0],str(con_end[0][1] + 1), str(con_start[0][1] + 1)]

				if line[5] in hg38_pos.keys():
					helper = hg38_pos[line[5]]
					helper.append([new_pos_l[0],new_pos_l[1], new_pos_l[2],  line[1], line[2], line[3], line[4]])
				else:
					hg38_pos[line[5]] = [[new_pos_l[0], new_pos_l[1], new_pos_l[2],  line[1], line[2],line[3], line[4]]] #chr start end 	major	minor	rsID	MAF	leadSNP_hg19
				o.write(new_pos + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\n')

	print("lost snps: " + str(counter_lost))
	return hg38_pos


## returns the numbe rof lines in a file without a header or a header starting with #
def getNumber(inputFile):

	counter = 0
	with open(inputFile) as i:
		for line in i:
			if line[0] != "#": # skip optional header of a bed file
				counter+= 1

	return counter


def generateRandomSNP(matchedSNPs):

	number = randrange(len(matchedSNPs))	
	info_  = matchedSNPs[number]
	rsID =  info_[5]
	if rsID == ".":
		rsID = "-"

	random_snp = info_[0] + '\t' + info_[1] + '\t' + info_[2] + '\t' + info_[3] + '\t' + info_[4] + '\t' + rsID + '\t' + info_[6] + '\n'

	return random_snp


if (len(sys.argv) < 2):
	print("python3 splitRandomSNPs.py, randomSNPDir")
else:

	randomSNPdir = sys.argv[1] ## result folder of the random snps from snpsnap and sneep
	# get the subfolders of the inputfolder
	subfolders = [ f.path for f in os.scandir(randomSNPdir) if f.is_dir() ]


	## lift SNPs to hg38 
	converter = get_lifter('hg19', 'hg38')

	## initalize a random seed for a reproducable result 
	random.seed(237)
	
	for elem in subfolders:
		celltype = elem.split("_")[-1]
		print(celltype)

		## lift position to hg38
		hg38_pos = lift_hg19_to_hg38(converter, elem + "/allRandomSNPs_withAlleleInfo_hg19.txt", elem + "/allRandomSNPs_withAlleleInfo_hg38.txt")

		## get number uniq snps original sneep input per celltype
		number_snps = getNumber("sneep_result_celltypeSpecific/sneep_cadlincGWAS_" + celltype + "_2.52x0.00001_noSampling/snpRegions.bed")

		seedSNPs = hg38_pos.keys()
		rounds = {} # number rounds -> snps
		
		## mk output dir 
		## store sampled snps in a extra directory but copy them afterwards to randomSNPs_sneep/ 
		if not os.path.isdir(os.path.join(elem, "randomSNPs_sneep/")):
			os.mkdir(elem + "/randomSNPs_sneep/")
		counter = 0
		# write for the 500 background rounds the matched snps seperated in an output file 
		for n in range(0,BACKGROUND_ROUNDS):

			uniq_snps = []
			## read in sneep sampled snps -> needs to be done first because we cannot resample them easily
			with open(elem +  "/sneep_cadlincGWAS_"+ celltype + "_2.52x0.00001_snpSampling_22_02/sampling/randomSNPs_" + str(n) + ".txt", 'r') as i, open(elem +  "/randomSNPs_sneep/randomSNPs_" + str(n) + ".txt", 'a') as o:
				for line in i:
					uniq_snps.append(line)

			## read in and check matched snps
			for s in seedSNPs: ## iterate over all snps and split them into 500 batches
				matchedSNPs = hg38_pos[s] ## get the matched snps
				if len(matchedSNPs) > n: ## write the top 500 to file 
					info_  = matchedSNPs[n]
					#write output in sneep format
					rsID =  info_[5]
					if rsID == ".":
						rsID = "-"
					current_snp = info_[0] + '\t' + info_[1] + '\t' + info_[2] + '\t' + info_[3] + '\t' + info_[4] + '\t' + rsID + '\t' + info_[6] + '\n'
					if not current_snp in uniq_snps:
						uniq_snps.append(current_snp)
					else:
						counter += 1		
						randomSNP = generateRandomSNP(matchedSNPs)
						while randomSNP in uniq_snps:
							randomSNP = generateRandomSNP(matchedSNPs)
						uniq_snps.append(randomSNP)
				else: ## what if there are not enough matched SNPs -> get a random one (which is a doublicated one)
					counter += 1		
					randomSNP = generateRandomSNP(matchedSNPs)
					while randomSNP in uniq_snps:
						randomSNP = generateRandomSNP(matchedSNPs)
					uniq_snps.append(randomSNP)

			if number_snps != len(uniq_snps): 
				print("ERROR " + celltype + " " + "round: " + str(n) + " " +  str(number_snps) + " " + str(len(uniq_snps)))
			with open(elem +  "/randomSNPs_sneep/randomSNPs_" + str(n)  + ".txt", 'w') as o:
				for u in uniq_snps:
					o.write(u)

		print("sampled SNP overall background rounds per celltype: " + str(counter))

