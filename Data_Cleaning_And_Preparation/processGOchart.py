#!/usr/bin/python
import sys

if len(sys.argv) == 1:
    print("processGOchart.py <GO Chart file>")
else:
    f = open(sys.argv[1],"r")
    data = f.read()
    f.close()
    data = data.split("\n")
    data = data[:len(data)-1]
    f = open(sys.argv[1].split(".")[0] + ".csv","w")
    f.write("GOClass,GOTerm,Count,P-value,Benjamini\n")
    for entry in data:
        terms = entry.split("\t")
        if terms[0][0:6] == "GOTERM":
            gotermentry = [terms[0],terms[1].replace(",","-"),terms[2],terms[4],terms[11]]
            f.write(",".join(gotermentry) + "\n")
    f.close()
