def main():
    fileIn1 = open("syndiniales_asvs.txt", "r")
    fileIn2 = open("dna-sequences.fasta", "r")
    fileOut = open("syndiniales_seqs.fasta","w")

    all = fileIn1.readlines()

    irofile = iter(fileIn2)
    for line in irofile:
        if ">" in line:
            seqName=line.replace(">","")
            for item in all:
                if item==seqName:
                    seqSeq = next(irofile)
                    print(">",seqName,sep="",end="",file=fileOut)
                    print(seqSeq,end="",file=fileOut)


    fileIn1.close()
    fileIn2.close()
    fileOut.close()




main()
