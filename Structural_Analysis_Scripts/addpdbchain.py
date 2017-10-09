#!usr/bin/python
import sys

if len(sys.argv) == 1:
    print "addpdbchain.py <amino acid sequence>"
else:
    output = ""
    f = open(sys.argv[1],"r")
    for line in f:
        if len(line) > 4:
            output += line[0:21] + "A" + line[22:]
        else:
            output += line
    f.close()
    f = open(sys.argv[1],"w")
    f.write(output)
    f.close()