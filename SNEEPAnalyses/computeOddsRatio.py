import sys, os
#from scipy.stats import binom_test
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

ROUNDS = 500 # number random background rounds snp
CUTOFF = 5 # odds-ratio cutoff
NUMBER_REGSNPs = 5
CUTOFF_P = 0.05 # pvalue cutoff loss and gain FDR corrected

if (len(sys.argv) < 6):
	print("python3 parseTFCount.py sneep_cell_type_dir, TFListFile, numberSNPsPerCelltypeFile, outputFile, outputFile2")
else:
	resultDir = sys.argv[1] #sneep_cell_types
	TFListFile = sys.argv[2] # list with the name of all tested TFs
	numberSNPsFile = sys.argv[3] # file with the number of SNPs per celltype numberSNPsPerCelltype.txt 
	outputFile = sys.argv[4] #heatmap
	outputFile2 =sys.argv[5] # summary file

	allTFs = []
	with open(TFListFile, 'r') as i:
		for line in i: 
			line = line.strip().split()
			allTFs.append(line[1])

    #gather number of snps per cell type
	celltypeSNPs = {} # celltype -> number SNP
	with open(numberSNPsFile, 'r') as i:
		for line in i: 
			line = line.strip().split('\t')
			helper = line[1].split("_")[5]
			celltypeSNPs[helper] = float(line[0])

    ## collect all celltypes and sort them
	celltypes = []
	for file_path in os.listdir(resultDir): ## for all celltypes
		celltype =  file_path.strip().split('_')[2]
		celltypes.append(celltype)

	celltypes.sort() # sort the celltypes alphabetically

	sigTFs = set()
	summary_oddsratio = {} #per celltype store oddsRatio
	summary_counts = {} # per cell type realCounts
	summaryGainLoss = {} # per cell type store gainLoss  pvalue FDR corrected
	summaryGainLoss_effect = {} # per cell type store gain loss effect
    ## compute odds ratio per celltyp

	consideredCelltypes = set()
	for c in celltypes:
		helper = {} # tf to realcount
		counts_random = {} # per TF the count over all background rounds
		odds_ratio = {} #per TF odds ratio
		with open(resultDir + "/sneep_cadlincGWAS_" + c + "_2.52x0.00001_backgroundSampling_22_02/TF_count.txt", 'r') as i:
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
            
            		#compute odds-ratio
			counter = 1
			numberSNPs = celltypeSNPs[c]
			for elem in realData:
				tf = header[counter]
				helper[tf] = float(elem)
				counter = counter + 1
				if float(elem) != 0:
					alpha = (float(elem) +1)/ numberSNPs # add pseudocount
					beta = counts_random[tf]/ numberSNPs
					odds_ratio[tf]  = (alpha / (1-alpha)) / (beta / (1-beta))
				else:		
					odds_ratio[tf]  = 0.0

		summary_oddsratio[c] = odds_ratio
		summary_counts[c] = helper

		## compute gain-loss info
		count_loss = {} # TFs: number of lost TF binding sites (negative DMax)
		count_gain = {} # TFs: number of gained TFBS (positive DMax)
		with open(resultDir + "/sneep_cadlincGWAS_" + c + "_2.52x0.00001_backgroundSampling_22_02/result.txt", 'r') as j:
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

		## compute for each TF two sided binomial test
		effect = 0
		pvalues = []
		gain_loss = {} # per TF is it a gain or a loss (gain = 1,loss = 0)
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
			if pvalue <= 0.05:
				print("celltype: " + c + " TF: " + elem + " countGain: " + str(gain) + " countLoss: " + str(loss) + " pvalue: " + str(pvalue) + " oddsRatio: " + str(odds_ratio[elem]))

			pvalues.append(pvalue)
		summaryGainLoss_effect[c] = gain_loss
        	##FDR correction
		test = multipletests(pvalues, alpha = 0.05, method = "fdr_bh", is_sorted = False, returnsorted = False) 
		corrected_pvalues = test[1] # grap only the adjusted pvalues
		gainLoss_corrected = {} # tf to corrected Pvalue
		counter = 0
		for elem in allTFs:
			gainLoss_corrected[elem] = corrected_pvalues[counter]
			counter += 1
		
		summaryGainLoss[c] = gainLoss_corrected	

		# get sig TFs over all celltypes
		counter = 0
		for tf in allTFs:
			if odds_ratio[tf] >= CUTOFF and  helper[tf] >= NUMBER_REGSNPs: ## found an interesting TF
				sigTFs.add(tf)
				consideredCelltypes.add(c)

	print(set(celltypes).difference(consideredCelltypes))
   	##write the header to the outputFiles
	with open(outputFile, 'w') as o1, open(outputFile2, 'w') as o2:
		o1.write("celltype") ## outputFile holds the sig odds-ratio for each TF over all celltypes
		o2.write("celltype") ## outputFile holding the corresponding gain/loss values
		#for tf in allTFs:
		for tf in sigTFs:
			o1.write("," + tf)
			o2.write("," + tf)
		o1.write("\n")
		o2.write("\n")

        ## write info to outputFile and store info for summaryFile
	with open(outputFile, 'a') as o1, open(outputFile2, 'a') as o2:
		for c in consideredCelltypes:
			o1.write(c)
			o2.write(c)
			counter = 0
			odds_ratio = summary_oddsratio[c]
			gainLoss = summaryGainLoss[c] 
			counts = summary_counts[c]
			effect = summaryGainLoss_effect[c]
			for tf in sigTFs:
				if odds_ratio[tf] >= CUTOFF and  counts[tf] >= NUMBER_REGSNPs: 
					o1.write("," + str(odds_ratio[tf])) ## write odds-ratio
					p = gainLoss[tf]
					if p <= CUTOFF_P:
						if effect[tf] == 1:
							o2.write(",1") ## it is a sig. gai
						else:
							o2.write(",-1") ## it is a sig. loss
					else:
						o2.write(",0") ## it is not sig.

				else:
					o1.write(",0.0") ## set odds-ratio to 2 
					o2.write(",2") ##  it is not relevant since TF is not sig associated to the celltype

				counter = counter + 1
			o1.write('\n')
			o2.write('\n')

	print("sigTFs: " + str(len(sigTFs)) + " " + str(sigTFs))
