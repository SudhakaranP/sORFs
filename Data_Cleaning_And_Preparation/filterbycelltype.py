#!/usr/bin/python
import sys

if len(sys.argv) == 1:
    print("filterbycelltype.py <classification> <fastafile> <genelist>")
else:
    f = open(sys.argv[2])
    fasta = f.read()
    f.close()
    f = open(sys.argv[3])
    genelist = f.read()
    f.close()

    fasta = fasta.split("\n")
    fasta = fasta[:len(fasta)-1]
    genelist = genelist.split("\n")
    genelist = genelist[:len(genelist)-1]

    genesubset = []

    for entry in fasta:
        if entry[0] == ">":
            temp = entry.split("|")
            if temp[9] == sys.argv[1]:
                for gene in genelist:
                    temp2 = gene.split(",")
                    if temp[4] == temp2[2]:
                        if (temp[5] < temp2[3] and temp[6] > temp2[3]) or (temp[5] > temp2[3] and temp[5] < temp2[4]):
                            if not(gene in genesubset):
                                genesubset.append(gene)

    f = open(sys.argv[3].split(".")[0] + sys.argv[1] + "subset.csv","w")
    for entry in genesubset:
        f.write(entry + "\n")
    f.close()