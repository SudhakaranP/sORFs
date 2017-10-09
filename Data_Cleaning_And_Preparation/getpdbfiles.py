#!/usr/bin/python
import sys
import requests

rcsbsearch = "https://www.rcsb.org/pdb/rest/search"

if len(sys.argv) == 1:
    print("getpdbfiles.py <pathtofile>")
else:
    output = ""
    f = open(sys.argv[1])
    data = f.read()
    f.close()
    data = data.split("\n")
    for l in data:
        xml = "<?xml version='1.0' encoding='UTF-8'?><orgPdbQuery><queryType>org.pdb.query.simple.UniprotGeneNameQuery</queryType><description>Uniprot Gene Name:"+ l + "</description>" +"<query>" + l + "</query>" + "</orgPdbQuery>"
        headers = {'Content-Type':'application/x-www-form-urlencoded'}
        result = requests.post(rcsbsearch,data=xml,headers=headers)
        if result.text != "":
            print(l)