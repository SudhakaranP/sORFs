#!/usr/bin/python
import sys
import os
import subprocess

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
TARGETP = "/home/david/Desktop/targetp-1.1/targetp"

if len(sys.argv) == 1:
    print("predictloc.py <data directory>")
else:
    f = open(WD + "/" + sys.argv[1] + "-loc.csv","w")
    f.write("Name,Len,mTP,SP,other,Loc,RC,TPlen\n")
    for root, dirs, files in os.walk(sys.argv[1]):
        for seqdirs in dirs:
            print("Processing " + seqdirs)
            args = [TARGETP, "-N", "-c", WD + "/" + sys.argv[1] + "/" + seqdirs + "/" + seqdirs + ".faa"]
            print(args)
            p = subprocess.check_output(args)
            p = p.split("\n")
            for line in p:
                if len(line) > 2:
                    if line[0:2] == "sO":
                        print(line)
                        output = ""
                        temp = line.split(" ")
                        for term in temp:
                            if len(term) > 0:
                                output += term + ","
                        output = output[:len(output)-1]
                        output += "\n"
                        f.write(output)
    f.close()
