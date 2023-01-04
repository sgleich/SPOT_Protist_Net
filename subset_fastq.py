def main():
    fileIn1 = open("try.txt", "r")
    fileIn2 = open("fastqtry.txt", "r")
    fileOut = open("syn_seqs.fasta","w")

    for i, line in enumerate(fileIn2):
        if i % 2 == 0:
            seq=line.replace(">","")
        if i % 2 ==1:
            seq2=line
            if seq in fileIn1:
                print(">",seq,sep="",end="",file=fileOut)
                print(seq2,end="",file=fileOut)


main()
