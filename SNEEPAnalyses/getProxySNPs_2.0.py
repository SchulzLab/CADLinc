import sys, os

def getProxy(output, chr_, snp_pos, snipaDir, r_2, current_snps):

	print(chr_)
	ld_file = snipaDir + "/grch37-1kgpp3v5-eur-chr" + str(chr_) + "-ld" #CHR	POS1(leadSNP pos)	POS2(proxySNP pos)	R2	D	DPRIME	RSID	RSALIAS	MINOR	MAF	MAJOR	CMMB	CM
	seen = set()
	with open(ld_file, 'r') as snipa: 
		snipa.readline() #skip header
		for line2 in snipa:
			line2 = line2.strip().split("\t")

			pos1 = line2[1]

			## a lead snp is found
			if pos1 in snp_pos:
				seen.add(pos1)
				proxyPos = line2[2]
				if pos1 == proxyPos:
					#write info from leadSNP file
					info = current_snps[pos1]
					for elem in info:
						o.write(elem + '\t' + line2[6] + '\t' + line2[9] + '\t' + "lead" +  '\t-\n') # add rsID and MAF, lead and lead_rsID
			
				else:
					current_r2 = float(line2[3])
					##real proxy snp found 
					if current_r2 >= r_2: 
						o.write("chr" + str(chr_) + "\t" + str(int(proxyPos) - 1)  + "\t" + proxyPos + '\t' +  line2[10] + "\t" + line2[8] + "\t" + line2[6] + "\t" + line2[9] + "\tproxy\t" + pos1 + "\n") # chr start end majorAllele minorAllel rsID MAF proxy lead_rsID

	return(snp_pos.difference(seen)) 

if (len(sys.argv) < 5):
	print("python3 getProxySNPs.py leadSNPsFile, r2, snipaDir, outputFile")
else:
	leadSNPsFile = sys.argv[1]
	r_2 = float(sys.argv[2])
	snipaDir = sys.argv[3]
	outputFile = sys.argv[4]

	chr_ = 1
	#chr_ = 21 ## for testing
	snp_pos  =  set()
	current_snps = {} #pos -> info from leadSNP file
	counter = 0
	with open(leadSNPsFile, 'r') as i, open(outputFile , 'w') as o:

		i.readline() #skip header
		for line in i: 

			line = line.strip().split('\t')	 #chr start end allele 1 allele2 maf pvalue
			pos = line[2] # start

			# still same chromosome?
			if line[0][3:] == str(chr_):
				# add pos to snp_pos
				snp_pos.add(pos)
				counter +=1
				if pos in current_snps.keys():
					helper = current_snps[pos]
					helper.append(line[0] + "\t" + line[1] + "\t" + line[2] + '\t' + line[3] + '\t' + line[4])
					#new da wird ja nichts angehangen
					current_snps[pos] = helper
				else:
					current_snps[pos] = [line[0] + "\t" + line[1] + "\t" + line[2] + '\t' + line[3] + '\t' + line[4]] # chr start end allele1 allele2
					
			## look up proxy snp for previous chromosome and add line as first entry in snp_pos and current_snps
			else:
				print("double check counter " + str(counter))
				print("number SNPs for chr" + str(chr_) + ": " + str(len(snp_pos)))
				notSeen = getProxy(o, chr_, snp_pos, snipaDir, r_2, current_snps)	
				print("snps not seen: " + str(len(notSeen)))

				for elem in notSeen:
					
					info = current_snps[elem]
					for j in info:
						o.write(j + '\t-\t-\tleadUnknown\t-\n') # add rsID and MAF  (lead Unknown and lead_rsID)

				#set new chr 
				chr_ += 1
				counter = 0
				snp_pos = set()
				snp_pos.add(pos)
				current_snps = {}
				current_snps[pos] = [line[0] + "\t" + line[1] + "\t" + line[2] + '\t' + line[3] + '\t' + line[4]]

		## for last chromosome 22
		print("number SNPs for chr" + str(chr_) + ": " + str(len(snp_pos)))
		print("double check counter " + str(counter))
		notSeen = getProxy(o, chr_, snp_pos, snipaDir, r_2, current_snps)	
		print("snps not seen " + str(len(notSeen)))

		for elem in notSeen:
					
			info = current_snps[elem]
			for j in info:
				o.write(j + '\t-\t-\tleadUnknown\t-\n') # add rsID and MAF

