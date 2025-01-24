import sys, os

if (len(sys.argv) < 3):
    print("python3 getUniqEntries_hg19.py swapped_leadAndProxy_hg19, outputFile")
else:
    inputSNPs = sys.argv[1]
    outputFile = sys.argv[2]

    considered = set() #chr:start:end:A1:A2 of valid SNPs
    not_considered = set()
    counter = 100
    mapping_status = {} # pos -> status  can be lead, proxy, leadAndAlsoProxy
    pos_rsID = {} # chr1:pos -> rsID

    c = 0
    ## write in file to get pos_rsID structure and to determine the status of a snp
    with open(inputSNPs, 'r') as i:
        for line in i:
            line = line.strip().split('\t') # chr1    2178998 2178999 C       T       rs4648818       0.38171 proxy   2164699
            a1 = line[3].upper()
            a2 = line[4].upper()
            snp = line[0] + ":" + line[1] + ":" + line[2] + ":" + a1 + ":" + a2 ## important to consider only uppercase alleles (otherwise file is not uniq)
            if (a1 == "A" or a1 == "C" or a1 == "G" or a1 == "T") and (a2 =="A" or a2 == "C" or a2 == "G" or a2 == "T"):
                considered.add(snp)
            else:
                if line[7] == "lead" or line[7] == "leadUnknown":
                    c +=1
            pos = line[0] + ":" + line[2]
            rsID = line[5]
            status = ""
            if snp in considered: 
                if not pos in pos_rsID.keys():
                    if rsID == "." or rsID == "-":  ## 10000 genome project knowns snp but does not have an rsID for it
                        rsID = "XX" + str(counter)
                        counter += 1
                    pos_rsID[pos] = rsID

                ## determine status
                if not snp in mapping_status.keys():
                    mapping_status[snp] = line[7] # set status to the one given in file
                else: # we have seen the snp before, is the status consistant
                    s = mapping_status[snp]
                    if s != line[7]:
                        mapping_status[snp] = "leadAndProxy"
            else:
                    not_considered.add(pos)
    print(c)
    print(" number considered snps: " + str(len(considered)))

    ## read in the file a second time to get the lead-proxy snp structure
    proxy_lead = {} # snps -> rsIDs of the correponding lead
    with open(inputSNPs, 'r') as i:
        for line in i:
            line = line.strip().split('\t') # chr1    2178998 2178999 C       T       rs4648818       0.38171 proxy   2164699
            snp = line[0] + ":" + line[1] + ":" + line[2] + ":" + line[3].upper() + ":" + line[4].upper() ## important to consider only uppercase alleles (otherwise file is not uniq)


            if line[8] != "-":
                pos_lead = line[0] + ":" + line[8]
                if not pos_lead in not_considered:
                    if snp in considered: 
                        if snp in proxy_lead.keys():
                            h = proxy_lead[snp]
                            h.add(pos_rsID[pos_lead]) ## add a new element to the set
                            proxy_lead[snp] = h 
                        else:
                            proxy_lead[snp] = {pos_rsID[pos_lead]} # define set
                            
                else:
                    if not snp in proxy_lead.keys():
                        proxy_lead[snp] = {"-"} # define set
            else:
                if not snp in proxy_lead.keys():
                    proxy_lead[snp] = {"-"} # define set

    seen_snps = set()
    with open(inputSNPs, 'r') as i, open(outputFile, 'w') as o:
        for line in i:
            line = line.strip().split('\t') # chr1    2178998 2178999 C       T       rs4648818       0.38171 proxy   2164699
            snp = line[0] + ":" + line[1] + ":" + line[2] + ":" + line[3].upper() + ":" + line[4].upper() ## important to consider only uppercase alleles (otherwise file is not uniq)
            pos = line[0] + ":" + line[2]
    
            if snp in considered and not snp in seen_snps: 
                seen_snps.add(snp)
                rsIDs = ""
                helper = proxy_lead[snp]
                for elem in helper:
                    rsIDs = rsIDs + elem + ","
                maf = line[6]
                if maf == "-" or maf == ".":
                    maf = -1 
                o.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3].upper() + '\t' + line[4].upper() + '\t' + pos_rsID[pos] +  '\t' + str(maf) + '\t' + rsIDs[:-1] + '\t'  +  mapping_status[snp] +  '\n')





            





            





