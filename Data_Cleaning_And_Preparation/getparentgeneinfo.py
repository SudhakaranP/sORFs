#!/usr/bin/python
import sys
import os
import requests
import xml.etree.ElementTree as ET

WD = "/media/david/genomicsdrive/NovelPeptideCharacterisation"
BASEURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

if len(sys.argv) == 1:
    print("getparentgeneinfo.py <data directory>")
else:
    for root, dirs, files in os.walk(sys.argv[1]):
        for seqdirs in dirs:
            os.chdir(WD + "/" + sys.argv[1] + "/" + seqdirs)
            print("Processing " + seqdirs)
            #Get parent gene id for current sequence
            f = open(seqdirs + ".fa")
            seqinfo = f.readline()
            f.close()
            seqinfo = seqinfo.split(" ")
            transcript_ids = seqinfo[3].split("=")[1].split(",")
            print(transcript_ids)
            #Access ncbi transcript records
            for transcript in transcript_ids:
                tparams = {"db" : "sequences", "id" : transcript, "retmode" : "xml"}
                transcript_request = requests.get(BASEURL, params=tparams)
                root = ET.fromstring(transcript_request.text)
                for qual in root.iter("GBQualifier"):
                    if qual.find("GBQualifier_name").text == "translation":
                        f = open(seqdirs + "-" + transcript + "-translation.faa","w")
                        f.write(">" + transcript + "-translation\n")
                        f.write(qual.find("GBQualifier_value").text)
                        f.close()

