#!/usr/bin/python
import sys
import os
import subprocess

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
MULTISCORE = "MultiLoc2-26-10-2009/src/multiloc2_prediction.py"

def readmsffile(f):
    msf = open(WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs +
                "-multiscore2.txt","r")
    data = msf.read()
    msf.close()
    data = data.split("\n")
    for line in data:
        if len(line) > 0:
            if line[0:4] == "gi":
                terms = line.split("\t")
                output = ""
                idterms = terms[0].split("|")
                output += idterms[1] + ","
                terms = terms[1:]
                outputlist = [0.0] * 9
                for term in terms:
                    locterms = term.split(":")
                    termloc = locterms[0]
                    outputlist[categories.index(termloc)] = locterms[1].lstrip()
                output += ",".join(outputlist) + "\n"
                print(output)
                f.write(output)

if len(sys.argv) == 1:
    print("runmultiscore2.py <data directory>")
else:
    categories = ["nuclear","cytoplasmic","peroxisomal","mitochondrial",
                    "Golgi apparatus","plasma membrane", "ER",
                    "extracellular","lysosomal"]
    seqcount = 0
    f = open(sys.argv[1] + "-multiscore.csv","w")
    f.write("ID," + ",".join(categories) + "\n")
    for root, dirs, files in os.walk(sys.argv[1]):
        for seqdirs in dirs:
            seqcount += 1
            print("Processing " + seqdirs + "---" + str(seqcount))
            if not os.path.isfile(WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs + "-multiscore2.txt"):
                args = ["python2", WD + "/" + MULTISCORE, "-fasta=" + WD + "/" +
                    sys.argv[1] + "/" + seqdirs + "/" + seqdirs + ".faa",
                    "-origin=animal", "-result=" + WD + "/" + sys.argv[1] +
                    "/" + seqdirs + "/" + seqdirs + "-multiscore2.txt"]
                p = subprocess.check_output(args)
            readmsffile(f)
    f.close()
