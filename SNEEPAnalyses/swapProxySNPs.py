import sys, os

if (len(sys.argv) < 4):
	print("python3 combineInfo.py summary_statistic, SNPFile, outputfile_hg19")
else:

	summaryStatisticFile = sys.argv[1]
	SNPFile = sys.argv[2]
	outputFileHg19 = sys.argv[3]

	leadSNPs = []
	proxySNPs = {} # pos -> [o_line, o_line]
	counter = 0
	with open(SNPFile, 'r') as i:
		for line in i:
			o_line = line
			line = line.strip().split('\t')
			if line[7] == "lead" or line[7] == "leadUnknown":
				counter = counter + 1
				leadSNPs.append(o_line)
	
			else:
				counter = counter + 1
				pos = line[0] + ":" +  line[2]
				if pos in proxySNPs.keys():
					h = proxySNPs[pos]
					h.append(o_line)	
					proxySNPs[pos] = h
				else:
					proxySNPs[pos] = [o_line]

	print(counter)

	seen_snps = []
	written_snps = set() # chr start end allele1 allele2 rsID
	with open(summaryStatisticFile, 'r') as i, open(outputFileHg19, 'w') as o:
		i.readline()
		for line in i:
			line = line.strip().split('\t')
			beta = float(line[9])
			cad_pos = "chr" + line[1] + ":" + line[2]
			a1 = line[3].upper()
			a2 = line[4].upper()
			if cad_pos in proxySNPs.keys():
				candidates = proxySNPs[cad_pos]

				for c in candidates:
					c_line = c.strip().split('\t')
					c_a1 = c_line[3].upper()
					c_a2 = c_line[4].upper()

					## check if alleles are the same
					if (c_a1 == a1 and c_a2 == a2) or (c_a1 == a2 and c_a2 == a1): ## alleles are the same
						if beta > 0.0:
							o.write(c_line[0] + '\t' + c_line[1] + '\t' + c_line[2] + '\t' + a1 + '\t' + a2 + '\t' + c_line[5] + '\t' + c_line[6] + '\t' + c_line[7] + '\t' + c_line[8] + '\n')
						else: # swap
							o.write(c_line[0] + '\t' + c_line[1] + '\t' + c_line[2] + '\t' + a2 + '\t' + a1 + '\t' + c_line[5] + '\t' + c_line[6] + '\t' + c_line[7] + '\t' + c_line[8] + '\n')
						## we were able two fix the allele positions of this snp
						## remember original snp
						written_snps.add(c_line[0]+ '\t' +  c_line[1] + '\t' +  c_line[2] + '\t' +  c_a1 + '\t' +  c_a2 + '\t' +  c_line[5])

		## write lead snps
		for elem in leadSNPs:
			o.write(elem)

		c_2 = 0	
		for elem in proxySNPs:
			helper = proxySNPs[elem]
			for h in helper: 
				h3 = h.split('\t')
				h2 = h3[0] + '\t' + h3[1] + '\t' + h3[2] + '\t' + h3[3] + '\t' + h3[4] + '\t' + h3[5]
				if not h2 in written_snps:
					c_2 = c_2 + 1
					o.write(h)

	print("counter proxy SNPs not found in summary statistic: " + str(c_2))
