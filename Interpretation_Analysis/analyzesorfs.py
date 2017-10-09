#!/usr/bin/python
import requests
import xml.etree.ElementTree as ET
import subprocess
import os

def pathexists(root, pathlist):
    for i in pathlist:
        if root.find(i) == None:
            return False
        else:
            root = root.find(i)
    return True

ncbigene = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene"
ncbinuccore = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
ncbiprotein = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"

f = open("data/PrabakaransORFsExonic.faa","r")
sorfs = f.read()
f.close()
sorfs = sorfs.split("\n")

f = open("exonicgenelistncbiid.csv","r")
exonicgenelist = f.read()
f.close()
exonicgenelist = exonicgenelist.split("\n")

genesorfmap = [] # ncbiid, gene name, sorfid, transcript id, cdsstart, cdsend, cdsframe, exons, sorfaliframe, sorfframe with reference to parent seq, rnaalistart, rnaaliend, rnaalieval, proteinid, proteinalistart, proteinaliend, proteinalieval

latestgene = ""
if os.path.isfile("PrabakaransORFsExonicParentGeneInfo.csv"):
    f = open("PrabakaransORFsExonicParentGeneInfo.csv","r")
    data = f.read()
    data = data.split("\n")
    data = data[:len(data)-1]
    f.close()
    latestgene = data[len(data)-1].split(",")[0]
    i = 0
    while i < len(data):
        if data[i].split(",")[0] == latestgene:
            del data[i]
        else:
            i += 1
    f = open("PrabakaransORFsExonicParentGeneInfo.csv","w")
    for entry in data:
        f.write(entry + "\n")
    f.close()
    outputf = open("PrabakaransORFsExonicParentGeneInfo.csv","a")
else:
    outputf = open("PrabakaransORFsExonicParentGeneInfo.csv","w")

for gene in exonicgenelist:
    if len(gene) > 0:
        if latestgene != "":
            if latestgene != gene.split(",")[0]:
                continue
            else:
                latestgene = ""
        temp = gene.split(",")
        for sorf in sorfs:
            if len(sorf) > 0:
                temp2 = sorf.split("|")
                if sorf[0] == ">":
                    if temp2[4] == temp[2]:
                        if (temp2[5] < temp[3] and temp2[6] > temp[3]) or (temp2[5] > temp[3] and temp2[5] < temp[4]):
                            print(temp[0])
                            result = requests.get(ncbigene + "&id=" + temp[0] + "&retmode=xml")
                            root = ET.fromstring(result.text)
                            geneproducts = []
                            products = {}
                            for child in root.find("Entrezgene").find("Entrezgene_locus"):
                                for geneproduct in child.find("Gene-commentary_products").findall("Gene-commentary"):
                                    geneproducts.append(geneproduct)
                            for product in geneproducts:
                                if product.find("Gene-commentary_type").text == "3":
                                    if pathexists(product, ["Gene-commentary_products", "Gene-commentary", "Gene-commentary_accession"]):
                                        products[product.find("Gene-commentary_accession").text] = product.find("Gene-commentary_products").find("Gene-commentary").find("Gene-commentary_accession").text
                                    else:
                                        products[product.find("Gene-commentary_accession").text] = ""
                                elif product.find("Gene-commentary_type").text == "22":
                                    products[product.find("Gene-commentary_accession").text] = ""
                            #find rna alignments
                            for rna, protein in products.items():
                                genesorfmapentry = [temp[0],temp[1], temp2[1], rna]
                                #get cds and exons
                                transcriptinfo = requests.get(ncbinuccore + "&id=" + rna + "&retmode=xml")
                                root = ET.fromstring(transcriptinfo.text)
                                exonlist = ""
                                cdsexist = False
                                for child in root.find("GBSeq").find("GBSeq_feature-table").findall("GBFeature"):
                                    if child.find("GBFeature_key").text == "CDS":
                                        cdsexist = True
                                        genesorfmapentry.append(child.find("GBFeature_intervals").find("GBInterval").find("GBInterval_from").text)
                                        genesorfmapentry.append(child.find("GBFeature_intervals").find("GBInterval").find("GBInterval_to").text)
                                        genesorfmapentry.append(str(int(genesorfmapentry[4]) % 3))
                                    elif child.find("GBFeature_key").text == "exon":
                                        exonlist += str(child.find("GBFeature_intervals").find("GBInterval").find("GBInterval_from").text)
                                        exonlist += "-"
                                        exonlist += str(child.find("GBFeature_intervals").find("GBInterval").find("GBInterval_to").text)
                                        exonlist += "|"
                                if not cdsexist:
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                if len(exonlist) > 0 :
                                    exonlist = exonlist[:len(exonlist)-1]
                                else:
                                    exonlist = "NA"
                                genesorfmapentry.append(exonlist)
                                #get transcript fasta
                                rnafasta = requests.get(ncbinuccore + "&id=" + rna + "&rettype=fasta")
                                f = open("PrabakaransORFsExonic/" + temp2[1] + "/" + rna + ".fna", "w")
                                f.write(rnafasta.text)
                                f.close()
                                #align transcript and sorf
                                significanthit = False
                                args = ["nhmmer", "--tblout", "PrabakaransORFsExonic/" + temp2[1] + "/" + rna +"-" + temp2[1] + ".txt",
                                    "PrabakaransORFsExonic/" + temp2[1] + "/" + temp2[1] + ".fna",
                                    "PrabakaransORFsExonic/" + temp2[1] + "/" + rna + ".fna"]
                                p = subprocess.check_output(args)
                                f = open("PrabakaransORFsExonic/" + temp2[1] + "/" + rna +"-" + temp2[1] + ".txt")
                                alignment = f.read()
                                f.close()
                                alignment = alignment.split("\n")
                                for line in alignment:
                                    if len(line) > 0:
                                        if line[0] != "#":
                                            significanthit = True
                                            aliparams = line.split(" ")
                                            aliparams = list(filter(None,aliparams))
                                            genesorfmapentry.append(str(int(aliparams[4]) % 3))
                                            genesorfmapentry.append(str(int(aliparams[6]) % 3))
                                            genesorfmapentry.append(aliparams[6])
                                            genesorfmapentry.append(aliparams[7])
                                            genesorfmapentry.append(aliparams[12])
                                            break
                                if not significanthit:
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                if protein != "":
                                    genesorfmapentry.append(protein)
                                    #get protein fasta
                                    proteininfo = requests.get(ncbiprotein + "&id=" + protein + "&rettype=fasta")
                                    f = open("PrabakaransORFsExonic/" + temp2[1] + "/" + protein + ".faa", "w")
                                    f.write(proteininfo.text)
                                    f.close()
                                    args = ["phmmer", "--domtblout", "PrabakaransORFsExonic/" + temp2[1] + "/" + protein +"-" + temp2[1] + ".txt",
                                    "PrabakaransORFsExonic/" + temp2[1] + "/" + temp2[1] + ".faa",
                                    "PrabakaransORFsExonic/" + temp2[1] + "/" + protein + ".faa"]
                                    p = subprocess.check_output(args)
                                    f = open("PrabakaransORFsExonic/" + temp2[1] + "/" + protein +"-" + temp2[1] + ".txt")
                                    alignment = f.read()
                                    f.close()
                                    alignment = alignment.split("\n")
                                    significanthit = False
                                    for line in alignment:
                                        if len(line) > 0:
                                            if line[0] != "#":
                                                significanthit = True
                                                aliparams = line.split(" ")
                                                aliparams = list(filter(None,aliparams))
                                                genesorfmapentry.append(aliparams[17])
                                                genesorfmapentry.append(aliparams[18])
                                                genesorfmapentry.append(aliparams[12])
                                                break
                                    if not significanthit:
                                        genesorfmapentry.append("NA")
                                        genesorfmapentry.append("NA")
                                        genesorfmapentry.append("NA")
                                else:
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                    genesorfmapentry.append("NA")
                                print(genesorfmapentry)
                                outputf.write(",".join(genesorfmapentry) + "\n")
outputf.close()





