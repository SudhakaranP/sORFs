#!/usr/bin/python
import sys
import requests
import xml.etree.ElementTree as ET

ncbigenesearchurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="

if len(sys.argv) == 1:
    print("sqltabletocsv.py <pathtofile>")
else:
    f = open(sys.argv[1])
    data = f.read()
    f.close()
    data = data.split("\n")
    data = data[:len(data)-1]
    genelist = []
    for gene in data:
        temp = gene.split("|")
        if len(genelist) > 0 and len(temp) == 6:
            latest = len(genelist)-1
            if genelist[latest][1] == temp[1].strip():
                if genelist[latest][3] > temp[3].strip():
                    genelist[latest][3] = temp[3].strip()
                if genelist[latest][4] < temp[4].strip():
                    genelist[latest][4] = temp[4].strip()
            else:
                result = requests.get(ncbigenesearchurl + temp[1].strip() + "[Gene Name]+AND+mus musculus[Organism]&retmode=xml")
                root = ET.fromstring(result.text)
                ncbiid = ""
                for child in root:
                    if child.tag == "IdList":
                        if len(child) > 0:
                            ncbiid += child[0].text
                generecord = [ncbiid,temp[1].strip(),temp[2].strip(),temp[3].strip(),temp[4].strip()]
                genelist.append(generecord)
                print(generecord)
        else:
            result = requests.get(ncbigenesearchurl + temp[1].strip() + "[Gene Name]+AND+mus musculus[Organism]&retmode=xml")
            root = ET.fromstring(result.text)
            ncbiid = ""
            for child in root:
                if child.tag == "IdList":
                    if len(child) > 0:
                        ncbiid += child[0].text
            generecord = [ncbiid,temp[1].strip(),temp[2].strip(),temp[3].strip(),temp[4].strip()]
            genelist.append(generecord)
            print(generecord)
    f = open(sys.argv[1].split(".")[0] + "ncbiid.csv","w")
    for generecord in genelist:
        f.write(generecord[0] + "," + generecord[1] + "," + generecord[2] + "," + generecord[3] + "," + generecord[4] + "\n")
    f.close()
