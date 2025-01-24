import sys, os
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

ROUNDS = 500 # number random background rounds snp
CUTOFF = 2 # odds-ratio cutoff
NUMBER_REGSNPs = 5 # number of reg SNPs per TF (try to avoid to inlcude TFs with low counts with high odds-ratio)
CUTOFF_P = 0.05 # pvalue cutoff loss and gain FDR corrected

if (len(sys.argv) < 6):
	print("python3 computeOddsRatio_mergedEnhancer.py sneep_cell_type_dir, TFListFile, numberSNPsPerCelltypeFile, outputFile_summryTFs, outputFile_loosGain")
else:
	resultDir = sys.argv[1] #sneep_cell_types
	TFListFile = sys.argv[2] # list with the name of all tested TFs
	numberSNPs = int(sys.argv[3]) # give the number of snps
	print(numberSNPs)
	outputFile_TFsummary = sys.argv[4]
	outputFile_lossGain = sys.argv[5]

	allTFs = []
	with open(TFListFile, 'r') as i:
		for line in i: 
			line = line.strip().split()
			allTFs.append(line[1])
	sigTFs = set()

	## collect infos for TF summary table
	TF_regSNPs = {} # list per TF the reg snps
	TF_genes = {} # list per TF the associated genes
	rsID_info = {} # rsID to chr:start-end_var1_var2
	mapping_geneNames = {} #ensemblID -> geneName
	with open(resultDir + "/result.txt", 'r') as j:
		j.readline()
		for r in j:
			r = r.strip().split('\t')
			rsID = r[3]
			## add rsID
			rsID_info[rsID] = r[0] + "_" + r[1] + "_" + r[2] #chr:start-end_A_C
			tf = r[6]
			## collect rsIDs over all SNPs
			if tf in TF_regSNPs.keys(): ## add a new rsID to a TF
				h = TF_regSNPs[tf]
				h.append(rsID)
				TF_regSNPs[tf] = h
			else:
				TF_regSNPs[tf] = [rsID]
			##collect genes over all SNPs 
			genes = r[15].strip().split(",")
			gene_names = r[16].strip().split(",")
			counter = 0
			for elem in genes: 
				mapping_geneNames[elem] = gene_names[counter]
				counter+=1 

			## add all genes of a TF
			if tf in TF_genes.keys(): 
				h = TF_genes[tf]
				for elem in genes:
					h.append(elem)
				TF_genes[tf] = h

			else:
				h = []
				for elem in genes:
					h.append(elem)
				TF_genes[tf] = h

	print(mapping_geneNames)
    	## compute odds ratio per celltyp
	counts_random = {} # per TF the count over all background rounds
	odds_ratio = {} #per TF odds ratio
	with open(resultDir + "/TF_count.txt", 'r') as i, open(outputFile_TFsummary, 'w') as o_1:
		header = i.readline() # first line is header
		header = header.strip().split('\t')

		realData = i.readline() # counts of the analysed SNPs
		realData = realData.strip().split('\t')[1:] # remove first entry from line 
		for line in i: # read in random data for the sampled rounds (500)
			counter = 1
			line = line.strip().split('\t')
			line = line[1:] # remove first entry which is the number of rounds
			for elem in line: 
				tf = header[counter]
				counter = counter + 1
				if tf in counts_random.keys():
					counts_random[tf] = counts_random[tf] + float(elem)
				else: 
					counts_random[tf] = float(elem)
		for k in counts_random.keys():
			counts_random[k] = counts_random[k]/ROUNDS # compute average random background signal 
			if counts_random[k] == 0.0:		
				counts_random[k] = 1

		counter = 1
		o_1.write("TF\tcount\tregSNPs_rsID\tregSNPs_info\toddsRatio\ttargetGenes\ttargetGenesName\texceedingCutoffs\n")
		for elem in realData:
			#print(elem)
			tf = header[counter]
			counter = counter + 1
			considered = False
			if float(elem) != 0:
				alpha = (float(elem) + 1) / numberSNPs ## add pseudocount
				beta = counts_random[tf]/ numberSNPs
				odds_ratio[tf]  = (alpha / (1-alpha)) / (beta / (1-beta))
				## get the TFs where the TF odds-ratio is higher 5 and more than 10 regSNPs are found per TF
				if odds_ratio[tf] >= CUTOFF and float(elem) >= NUMBER_REGSNPs:
					considered = True
					sigTFs.add(tf)
				## summaries rsIDs and the snp info
				s = "" 
				i_ = ""
				seen = set()
				for rSNP in TF_regSNPs[tf]:
					if not rSNP in seen:
						s += rSNP + ","
						i_ += rsID_info[rSNP] + ","
						seen.add(rSNP)

				## summaries genes and the gene name
				g = ""
				g_n = ""
				for gene in TF_genes[tf]:
					if not gene in seen:
						g +=gene + ","
						g_n += mapping_geneNames[gene] + ","
						seen.add(gene)
	
				o_1.write(tf + '\t' + elem + '\t' + s[:-1] + '\t' + i_[:-1] + '\t' +   str(odds_ratio[tf]) + '\t' + g[:-1] + '\t' + g_n[:-1] + '\t'  + str(considered) + "\n")
			else:		
				odds_ratio[tf]  = 0.0
				o_1.write(tf + '\t' + elem + '\t' + "-" + '\t' + "-" + '\t' +   str(odds_ratio[tf]) + '\t' + "-" + '\t' + "-" + '\t'  + str(considered) + "\n")

	## compute gain-loss info
	count_loss = {} # TFs: number of lost TF binding sites (negative DMax)
	count_gain = {} # TFs: number of gained TFBS (positive DMax)
	with open(resultDir + "/result.txt", 'r') as j:
		j.readline()
		for r in j:
			r = r.strip().split('\t')
			tf = r[6]
			score = float(r[12])
			if score < 0.0: # loss binding site
				if tf in count_loss.keys():
					count_loss[tf] =  count_loss[tf] + 1 
				else:
					count_loss[tf] = 1
			if score > 0.0: # gain binding site 
					if tf in count_gain.keys():
						count_gain[tf] =  count_gain[tf] + 1 
					else:
						count_gain[tf] = 1
	effect = 0
	pvalues = []
	gain_loss = {} # per TF is it a gain or a loss (gain = 1,loss = 0)
	with open(outputFile_lossGain, 'w') as o_1:
		o_1.write("TF\tcountGain\tcountLoss\tpvalue\toddsRatio\n") 
		for elem in allTFs: 
			gain = 0
			loss = 0
			pvalue = 1
			if elem in count_loss.keys():
				loss = count_loss[elem]
				if elem in count_gain.keys(): # otherwise the count for gains is 0
					gain = count_gain[elem]
				result = binomtest(k = gain, n = gain + loss, p =0.5, alternative = "two-sided")
				pvalue = result.pvalue
				if gain < loss: # loss
					gain_loss[elem] = 0
				else: # gain
					gain_loss[elem] = 1
			elif elem in count_gain.keys():
				gain = count_gain[elem]
				if elem in count_loss.keys(): # otherwise the count for gains is 0
					loss = count_loss[elem]
				result = binomtest(k = gain, n = gain + loss, p =0.5, alternative = "two-sided")
				pvalue = result.pvalue
				if gain < loss: # loss
					gain_loss[elem] = 0
				else: # gain
					gain_loss[elem] = 1
			else:
				pvalue = 1
				if gain < loss: # loss
					gain_loss[elem] = 0
				else: # gain
					gain_loss[elem] = 1
			o_1.write(elem + "\t" + str(gain) + "\t" + str(loss) + "\t" + str(pvalue) + "\t" + str(odds_ratio[elem]) + '\n')

			pvalues.append(float(pvalue))
        	##FDR correction
	test = multipletests(pvalues, alpha = 0.05, method = "fdr_bh", is_sorted = False, returnsorted = False) 
	corrected_pvalues = test[1] # grap only the adjusted pvalues

