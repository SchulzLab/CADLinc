import sys, os

if (len(sys.argv) < 5):
	print("python3 parseSNPSNAP_output.py ")
else:
	GWAS_snps_hg38 = sys.argv[1] # GWAS_summary_filtered_betaSwapped_2.52x0.00001_hg38_16_02_lead_proxySNPs_swapped_uniq.txt
	sneep_result_dir = sys.argv[2]
	matchedSNPsFile = sys.argv[3] # snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_16_02_matched_snps_for_input_snps/matched_snps.tsv
	qualityFile = sys.argv[4] # snpsnap_results/GWAS_cadlinc_2.52x0.00001_hg19_16_02_matched_snps_for_input_snps/match_quality_per_snp.tsv
	outputDir = sys.argv[5]

	## read in mapping between hg19 and hg38  
	hg19_hg38 = {} ## chr_pos (1based) hg19 to  chr_start_end hg38 (0-based)
	hg38_hg19 = {} ## chr_start_end hg38 (0-based) to chr_pos (1based)
	hg38_info = {} ## chr_start_end hg38 (0-based) to all info of this file

	with open(GWAS_snps_hg38) as i:
		i.readline() #skip header
		for line in i:
			line = line.strip().split('\t')
			pos_hg38 = line[0] + ":" + line[1] + "-" + line[2] + "_" + line[3] + "_" + line[4]
			helper = line[9].split(":") # split position such that we get 1:pos (1-based)
			pos_hg19 = helper[0][3:] + ":" + helper[1].split("-")[1]
			hg19_hg38[pos_hg19] = pos_hg38
			hg38_hg19[pos_hg38] = pos_hg19
			hg38_info[pos_hg38] = line

	print("read in mapping")
	## read in matched SNPs 
	matchedSNPs = {} # hg19 pos snp to list of matched SNPs
	with open (matchedSNPsFile, 'r') as i:
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			matchedSNPs[line[0]] = line[1:]
	print("read in matched SNPs")
			
	## read in quailty per snp
	qualitySNPs = {} # pos hg19 to num_uniqmatches
	with open(qualityFile, 'r') as i:
		i.readline() # skip header
		for line in i:
			line = line.strip().split('\t')
			qualitySNPs[line[0]] = int(line[1])
	print("read in quality SNPs")
	
	numberSNPs = set()	
	with open(outputDir + "summary.txt", 'w') as o:
		o.write("numberSNPsSneepInput\tnumberSNPsToLessMatches\tnumberSNPsNotMatched\n")
		
		## get for the current cell type the considered SNPs
		## when rems are used
		## when no enhancer-geen interactios are used and no footprint
		with open(sneep_result_dir + "/SNPsUnique.bed", 'r') as i:
			for line in i:
				line = line.strip().split('\t')
				## when no enhancer-geen interactios are used and no footprint
				pos_hg38 = line[0] + ":" + line[1] + "-" + line[2] + "_" + line[3]  + "_" + line[4]
				numberSNPs.add(pos_hg38)
		## build random SNP set per cell type
		matchedSNPs_ = [] #hg19
		c_1 = 0
		c_2 = 0
		if not os.path.isdir(os.path.join(outputDir, "randomSNPs")):
			os.mkdir(os.path.join(outputDir, "randomSNPs")) 
		with open(outputDir + "/randomSNPs/notMatchedSNPs.txt", 'w') as o_c:
			for snp in numberSNPs: # snp position is hg38
				pos_hg19 = hg38_hg19[snp] # get corresponding hg19 position
				if pos_hg19 in qualitySNPs.keys(): # the snp was found in their data
					matches = qualitySNPs[pos_hg19]
					if matches >50:
						matchedSNPs_.append(pos_hg19)
					else:
						info = hg38_info[snp]
						rsID = info[5]
						if rsID == "-":
							rsID = "."
						MAF = info[6]
						if MAF == "-":
							MAF = -1
						o_c.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + info[3] + '\t' + info[4] + '\t' + rsID + '\t' + str(MAF) + "\ttoLess\n")
						c_1 += 1

				else:
					info = hg38_info[snp]
					rsID = info[5]
					if rsID == "-":
						rsID = "."
					MAF = info[6]
					if MAF == "-":
						MAF = -1
					o_c.write(info[0] + '\t' + info[1] + '\t' + info[2] + '\t' + info[3] + '\t' + info[4] + '\t' + rsID + '\t' + str(MAF) + "\tnotListed\n")
					c_2 += 1


		with open(outputDir + "/randomSNPs/allRandomSNPsPos_hg19.txt", 'w') as o3:
			for elem in matchedSNPs_: #hg19
				r_snps = matchedSNPs[elem]
				for s in r_snps:
					o3.write(elem + '\t' + s + '\n')

		o.write(str(len(numberSNPs)) + '\t' + str(c_1) + '\t' + str(c_2) + '\n' )
	

