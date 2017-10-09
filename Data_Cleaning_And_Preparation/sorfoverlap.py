#!/usr/bin/python
import sys

if len(sys.argv) != 5:
    print("sorfoverlap.py <sorfs file> <chromosome> <start> <end>")
else:
    start = sys.argv[3]
    end = sys.argv[4]
    chromosome = sys.argv[2]
    f = open(sys.argv[1],"r")
    sorfs = f.read()
    f.close()
    output = ""
    addseq = False
    sorfs = sorfs.split("\n")
    for sorf in sorfs:
        if len(sorf) > 0:
            if sorf[0] == ">":
                temp = sorf.split("|")
                if chromosome == temp[4]:
                    if temp[5] >= start and temp[5] < end and temp[6] > start and temp[6] < end:
                        addseq = True
                        output += sorf + "\n"
            elif addseq:
                addseq = False
                output += sorf + "\n"
    f = open("overlappingsorfs.faa","w")
    f.write(output)
    f.close()

