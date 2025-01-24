import sys, os


def getProxy(snv_structure, chr_, snp_pos, snipaDir, r_2, mapping_rsID):

	print(chr_)
	ld_file = snipaDir + "/grch37-1kgpp3v5-eur-chr" + str(chr_) + "-ld" #CHR	POS1(leadSNP pos)	POS2(proxySNP pos)	R2	D	DPRIME	RSID	RSALIAS	MINOR	MAF	MAJOR	CMMB	CM
	seen = set()

	counter = 0
	with open(ld_file, 'r') as snipa: 
		snipa.readline() #skip header
		for line2 in snipa:
			line2 = line2.strip().split("\t")

			pos1 = line2[1]
			current_r2 = float(line2[3])
			if current_r2 >= r_2:  ## if r2 is too small it does not matter whether it is a proxy SNV or not 
				## a lead snp is found
				if pos1 in snp_pos:
					proxyPos = line2[2]
					if pos1 == proxyPos: ## SNP is to itself of course in LD -> but this information we do not  need
						continue	
					else: ##real proxy snp found 
						if proxyPos in snp_pos: ## is this proxy SNV in our SNV set, if so keep info else skip the SNV
							seen.add(pos1)
							counter += 1
							pos_0based =  "chr" + str(chr_) + ":" + str(int(pos1)-1) + "-" + str(pos1) ## get position chr:start-end og the current lead SNV in hg19
							if pos_0based in snv_structure:

								helper = snv_structure[pos_0based]["proxySNPs"] # contains a rsID: rsID, proxySNPs = {chr:start-end}
								helper.append(mapping_rsID[proxyPos])
								snv_structure[pos_0based]["proxySNPs"]  = helper
							else:
								snv_structure[pos_0based]  = {"proxySNPs" : [mapping_rsID[proxyPos]]}
	print("found LD SNVs: " + str(counter)) 

	return(snp_pos.difference(seen)) 

if (len(sys.argv) < 5):
	print("python3 getProxySNPs.py leadSNPsFile, r2, snipaDir, outputFile")
else:
	leadSNPsFile = sys.argv[1]
	r_2 = float(sys.argv[2])
	snipaDir = sys.argv[3]
	outputFile = sys.argv[4]

	chr_ = 1
	snp_pos  =  set()
	counter = 0
	snv_structure = {} # holds for all SNVs the SNVs in LD
	mapping = {} ## hg19 to HG38
	mapping_rsID = {} ## gives for each position the rsID
	with open(leadSNPsFile, 'r') as i:
		i.readline() #skip header

		for line in i: 
			o_line = line ## keep original line
			line = line.strip().split('\t')	 ##chr	start	end	allele1	allele2	rsId	MAF	rsIDsLead	status	position_hg19 (chr:start-end)
			
			#get hg19 pos 
			pos =  line[9].split("-")[1] #1-based
			start = line[9].split("-")[0].split(":")[1]

			## mapping hg19 position to coresponding lines
			if line[9] in mapping.keys():
				helper = mapping[line[9]] 
				helper.append(o_line)
				mapping[line[9]]  = helper
			else:
				mapping[line[9]]  = [o_line]

			## add position to mapping_rsID

			# still same chromosome?
			if line[0][3:] == str(chr_):

				mapping_rsID[pos] = line[5]
				if (abs(int(start) - int(pos)) == 1):
					# add pos to snp_pos
					snp_pos.add(pos)
					counter +=1
				else:
					print(line)
					
			## look up proxy snp for previous chromosome and add line as first entry in snp_pos and current_snps
			else:
				print("double check counter " + str(counter))
				print("number SNPs for chr" + str(chr_) + ": " + str(len(snp_pos)))
				notSeen = set()
				if len(snp_pos) != 0:
					notSeen = getProxy(snv_structure, chr_, snp_pos, snipaDir, r_2, mapping_rsID)
					print("snps not seen: " + str(len(notSeen)))
	
				for elem in notSeen:
					pos_hg19 = "chr" + str(chr_) + ":" +  str(int(elem) - 1) +  "-"  + str(elem)
					snv_structure[pos_hg19] = {"proxySNPs" : []}

				chr_ += 1

				counter = 0
				snp_pos = set()
				snp_pos.add(pos)
				mapping_rsID = {}
				mapping_rsID[pos] = line[5]

		## for last chromosome 22
		print("number SNPs for chr" + str(chr_) + ": " + str(len(snp_pos)))
		print("double check counter " + str(counter))
		notSeen = getProxy(snv_structure, chr_, snp_pos, snipaDir, r_2, mapping_rsID)	
		print("snps not seen " + str(len(notSeen)))

		for elem in notSeen:
			pos_hg19 = "chr" + str(chr_) + ":" +  str(int(elem) - 1) +  "-"  + str(elem)
			snv_structure[pos_hg19] = {"proxySNPs" : []}

	##write output in same format as before
	print(len(mapping.keys()))
	with open(outputFile, 'w') as o:
		o.write("chr\tstart\tend\tallele1\tallele2\trsId\tMAF\tldSNVs\tstatus\tposition_hg19\n")
		for elem in mapping.keys(): ## iterate over all SNVs
			ld_snvs = snv_structure[elem]["proxySNPs"] ##get LD SNVs
			ld_snvs_string = ""
			for s in ld_snvs:
				ld_snvs_string += s + ","

			ld_snvs_string = ld_snvs_string[:-1] # remove last ,
			
			snvs = mapping[elem]
			for line in snvs:
				line = line.strip().split('\t') ## list that contains the lines of the oinput file that are linked to the hg19 position
				o.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + ld_snvs_string + '\t' + line[8] + '\t' +  line[9] + '\n') ##chr  start   end     allele1 allele2 rsId    MAF     rsIDsLead       status  position_hg19 (chr:start-end)
	



			

					
