#!/usr/bin/python
import sys

def readcsv(filepath):
    f = open(filepath,"r")
    coords = f.read()
    f.close()
    coords = coords.split("\n")
    coords = coords[:len(coords)-1]
    return coords

def checkintersect(start1, end1, start2, end2):
    if ((start1 < start2) and (end1 >= start2) or (start1 > start2 and start1 < end2)):
        return True
    else:
        return False

if len(sys.argv) != 3:
    print("findintersect.py <coords1.txt> <coords2.txt>")
else:
    coords1 = readcsv(sys.argv[1])
    coords2 = readcsv(sys.argv[2])
    output = open(sys.argv[1].split(".")[0] + "-" + sys.argv[2].split(".")[0] + "-intersect.csv","w")
    for entry in coords1:
        terms1 = entry.split(",")
        for entry2 in coords2:
            terms2 = entry2.split(",")
            if terms1[0] == terms2[0]:
                if checkintersect(int(terms1[1]),int(terms1[2]),int(terms2[1]),int(terms2[2])):
                    output.write(",".join(terms1) + "," + ",".join(terms2) + "\n")
    output.close()