#!/usr/bin/python
import sys
import os

PATH_TO_DATA = "/media/david/genomicsdrive/NovelPeptideCharacterisation/data"
PATH_TO_OUTPUT = "/media/david/genomicsdrive/NovelPeptideCharacterisation"

if len(sys.argv) == 1:
    print("Please include a fasta file")
else:
    curseq = ""
    count = 0
    if not os.path.isdir(PATH_TO_OUTPUT + "/" +sys.argv[1].split(".")[0]):
        os.mkdir(PATH_TO_OUTPUT + "/" + sys.argv[1].split(".")[0], 0755)
    PATH_TO_OUTPUT = PATH_TO_OUTPUT + "/" + sys.argv[1].split(".")[0]
    with open(PATH_TO_DATA + "/" + sys.argv[1], "r") as seqfile:
        header = ""
        for line in seqfile:
            if line[0] == ">":
                header = line
                curseq = line.split("|")[1]
                if not os.path.isdir(PATH_TO_OUTPUT + "/" + curseq):
                    os.makedirs(PATH_TO_OUTPUT + "/" + curseq)
                print curseq
            else:
                newseqfile = open(PATH_TO_OUTPUT + "/" + curseq + "/" + curseq
                + ".faa", "w")
                newseqfile.write(header)
                newseqfile.write(line)
                newseqfile.close()
                print "done"
                count = count + 1
    print(str(count) + " Sequences processed")


