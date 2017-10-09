#!/usr/bin/python
import sys
import os
import subprocess

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
INTERPRO = "interproscan-5.25-64.0/interproscan.sh"

if len(sys.argv) == 1:
    print("runmultiscore2.py <data directory>")
else:
    seqcount = 0
    for root, dirs, files in os.walk(sys.argv[1]):
        for seqdirs in dirs:
            seqcount += 1
            print("Processing " + seqdirs + "---" + str(seqcount))
            if not os.path.isfile(WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs + "-interpro.tsv"):
                args = [WD + "/" + INTERPRO, "-i", WD + "/" +
                    sys.argv[1] + "/" + seqdirs + "/" + seqdirs + ".faa",
                    "-f", "tsv", "-goterms",
                    "-o", WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs + "-interpro.tsv"]
                p = subprocess.check_output(args)
                print(p)