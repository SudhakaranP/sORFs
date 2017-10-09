#!/usr/bin/python
import sys

PATH_TO_DATA = "/media/david/genomicsdrive/NovelPeptideCharacterisation/data"

if len(sys.argv) == 1:
    print("Please include a fasta file")
else:
    f = open(sys.argv[1].split(".")[0] + "-seqlen.csv","w")
    f.write("ID,LENGTH\n")
    curseq = ""
    count = 0
    with open(PATH_TO_DATA + "/" + sys.argv[1], "r") as seqfile:
        header = ""
        for line in seqfile:
            if line[0] == ">":
                header = line
                curseq = line.split("|")[1]
                print curseq
            else:
                f.write(curseq + "," + str(len(line)-1) + "\n")
                print "done"
                count = count + 1
    print(str(count) + " Sequences processed")
    f.close()

