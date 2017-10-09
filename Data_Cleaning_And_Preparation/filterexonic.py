#!/usr/bin/python
import sys

if len(sys.argv) == 1:
    print("filterexonic.py <pathtofile>")
else:
    f = open(sys.argv[1])
    data = f.read()
    f.close()
    data = data.split("\n")
    output = "select distinct name2, chrom, txStart, txEnd from refGene where"
    for i in range(0,len(data)):
        if len(data[i]) > 1:
            if data[i][0] == ">":
                temp = data[i].split("|")
                outputstr = ' (chrom="' + temp[4] + '" AND ((txStart < ' + temp[5] + ' AND txEnd > ' + temp[5] + ') OR (txStart < ' + temp[6] + ' AND txStart >= ' + temp[5] + '))) OR'
                output += outputstr
    output = output[:len(output)-3] + ";"
    f = open("exoniccommand.txt","w")
    f.write(output)
    f.close()
