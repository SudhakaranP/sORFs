#!/usr/bin/python
import sys

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
PSSCAN = "ps_scan/ps_scan.pl"
PRORULES = "ps_scan/prorules.txt"
if len(sys.argv) == 1:
    print("countprosite.py <data directory>")
else:
    f = open(WD + "/ps_scan/prorules.txt","r")
    prorules = f.read()
    f.close()
    f = open(sys.argv[1] + "-prosite.csv","r")
    data = f.read()
    f.close()
    data = data.split("\n")
    data = data[1:]
    f = open(sys.argv[1] + "-prosite.csv","w")
    f.write("ID," + prorules + "\n")
    prorules = prorules.split(",")
    for line in data:
        terms = line.split(",")
        seqid = terms[0]
        print("Processing " + seqid)
        pfile = open(WD + "/" + sys.argv[1] + "/" +
                    "/" + seqid + "/" + seqid + "_prosite.txt","r")
        pfiledata = pfile.read()
        pfile.close()
        pfiledata = pfiledata.split("\n")
        pterms = ["0"] * len(prorules)
        count = 0
        curterm = ""
        textterms = []
        for pline in pfiledata:
            if len(pline) > 0:
                if pline[0] == ">":
                    if curterm != "":
                        pterms[prorules.index(curterm)] = str(count)
                        textterms.append(curterm)
                        curterm = pline.split(" ")[3]
                    else:
                        curterm = pline.split(" ")[3]
                    count = 0
                else:
                    count += 1
        if curterm != "":
            pterms[prorules.index(curterm)] = str(count)
            textterms.append(curterm)
        if len(textterms) > 0:
            print textterms
        f.write(seqid + "," + ",".join(pterms) + "\n")
    f.close()
