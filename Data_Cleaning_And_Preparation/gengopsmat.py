#!/usr/bin/python
import os

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation/ps_scan"

os.chdir(WD)
f = open("prositegoterms.txt")
goterms = f.read()
f.close()
mat = {}
goterms = goterms.split(",")
curterm = ""
with open("prorule.dat", "r") as f:
    for line in f:
        if line[0:12] == "TR   PROSITE":
            temp = line.split(";")
            curterm = temp[2].lstrip()
            print curterm
            if curterm not in mat:
                mat[curterm] = ["0"] * (len(goterms) + 1)

        elif line[0:2] == "GO":
            temp = line.split(";")
            go = temp[0][5:]
            for i in range(0,len(goterms)):
                if goterms[i] == go:
                    mat[curterm][i+1] = "1"
                    break
f.close()
f = open("gopspair.csv","w")
f.write("ID," + ",".join(goterms) + "\n")
for k,v in mat.iteritems():
    f.write(k + "," + ",".join(v) + "\n")
f.close()

