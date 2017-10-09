#!/usr/bin/python
import sys
import os
import subprocess

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
PSSCAN = "ps_scan/ps_scan.pl"
PRORULES = "ps_scan/prorules.txt"

def readprositefile(f):
    pfile = open(WD + "/" + sys.argv[1] + "/" +
                "/" + seqdirs + "/" + seqdirs + "_prosite.txt","r")
    p = pfile.read()
    pfile.close()
    print(p)
    p = p.split("\n")
    terms = ["0"] * len(prorules)
    textterms = []
    curterm = ""
    count = 0
    for line in p:
        if len(line) > 0:
            if line[0] == ">":
                if curterm != "":
                    terms[prorules.index(curterm)] = str(count)
                    textterms.append(curterm)
                curterm = line.split(" ")[3]
                count = 0
            else:
                count += 1
    if curterm != "":
        terms[prorules.index(curterm)] = str(count)
        textterms.append(curterm)
    if len(textterms) > 0:
        print(textterms)
    f.write(seqdirs + "," + ",".join(terms) + "\n")

if len(sys.argv) == 1:
    print("prositescan.py <data directory>")
else:
    f = open(WD + "/ps_scan/prorules.txt","r")
    prorules = f.read()
    f.close()
    f = open(sys.argv[1] + "-prosite.csv","w")
    f.write("ID," + prorules + "\n")
    prorules = prorules.split(",")
    seqcount = 0
    for root, dirs, files in os.walk(sys.argv[1]):
        for seqdirs in dirs:
            seqcount += 1
            print("Processing " + seqdirs + "---" + str(seqcount))
            if not os.path.isfile(WD + "/" + sys.argv[1] + "/" + "/" + seqdirs + "/" + seqdirs + "_prosite.txt"):
                args = ["perl", WD + "/" + PSSCAN, "-d", WD + "/ps_scan/prosite.dat",
                        WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs + ".faa"]
                p = subprocess.check_output(args)
                pfile = open(WD + "/" + sys.argv[1] + "/" +
                "/" + seqdirs + "/" + seqdirs + "_prosite.txt","w")
                pfile.write(p)
                pfile.close()
            readprositefile(f)
    f.close()

