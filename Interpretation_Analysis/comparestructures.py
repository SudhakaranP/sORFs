#!/usr/bin/python
import os
import sys
from pymol import *
from shutil import copyfile

def getTopStructure(filename):
    f = open(filename,"r")
    data = f.readline()
    data = f.readline()
    result = {}
    f.close()
    if filename.split(".")[1] == "csv":
        data = data.split(",")
        result["Filename"] = data[0].split("/")[len(data[0].split("/"))-1]
        result["Score"] = data[5]
    elif filename.split(".")[1] == "txt":
        data = data.split("\t")
        result["Filename"] = data[0]
        result["Score"] = data[10]
    return result

aliases = {}
aliases["sORF0000440275"] = "ncORF0000440277"
aliases["sORF0000440394"] = "ncORF0000440398"
aliases["sORF0000440911"] = "ncORF0000440918"
aliases["sORF0000441046"] = "ncORF0000441054"
aliases["sORF0000442433"] = "ncORF0000442466"
aliases["sORF0000442601"] = "ncORF0000442636"
aliases["sORF0000442649"] = "ncORF0000442684"
aliases["sORF0000442714"] = "ncORF0000442753"


f = open("PrabakaransORFs2" + "-" + "PrabakaransORFs2refseq" + "-compare.csv","w")
f.write("sorfid,struct1,alphabetascore1,database1,foldingsoftware1,struct2,alphabetascore2,database2,foldingsoftware2,rmsd,couplingsprob,couplingsrank,normalisedcouplingrank,skippedcouplings\n")
for root, dirs, files in os.walk("PrabakaransORFs2"):
    for seqdirs in dirs:
        if "ORF" not in seqdirs:
            break
        entry = seqdirs + ","
        if os.path.isfile(root + "/" + seqdirs + "/ncORF" + seqdirs[4:] + ".faa"):
            prefix = "ncORF" + seqdirs[4:]
        elif seqdirs in aliases:
            if os.path.isfile(root + "/" + seqdirs + "/" + aliases[seqdirs] + ".faa"):
                prefix = aliases[seqdirs]
        else:
            prefix = seqdirs
        print("Processing - " + prefix)
        struct1 = None
        struct2 = None
        struct1couplings = None
        struct2couplings = None
        struct1couplingsrank = None
        struct2couplingsrank = None
        if os.path.isfile(root + "/" + seqdirs + "/fold/" + prefix + "_ranking.csv"):
            structinfo = getTopStructure(root + "/" + seqdirs + "/fold/" + prefix + "_ranking.csv")
            entry += structinfo["Filename"] + "," + structinfo["Score"] + ",uniprot100,local,"
            struct1 = root + "/" + seqdirs + "/fold/" + structinfo["Filename"]
            copyfile(struct1,"/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2TopPDBEVfoldranking/" + structinfo["Filename"])
            cmd.load(struct1,"struct1")
            cmd.hide("lines","all")
            cmd.show("cartoon","all")
            cmd.color("green","struct1")
            cmd.orient("struct1")
            cmd.png("96sORFstructurecomparison/" + structinfo["Filename"] + ".png", dpi=300)
            cmd.delete("struct1")
        elif os.path.isfile(root + "/" + seqdirs + "/structure_outputs/" + prefix + "_alphabeta_ranking.txt"):
            structinfo = getTopStructure(root + "/" + seqdirs + "/structure_outputs/" + prefix + "_alphabeta_ranking.txt")
            entry += structinfo["Filename"] + "," + structinfo["Score"] + ",uniprot100,server,"
            struct1 = root + "/" + seqdirs + "/structure_outputs/" + structinfo["Filename"] + ".pdb"
            copyfile(struct1,"/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2TopPDBEVfoldranking/" + structinfo["Filename"] + ".pdb")
            cmd.load(struct1,"struct1")
            cmd.hide("lines","all")
            cmd.show("cartoon","all")
            cmd.color("green","struct1")
            cmd.orient("struct1")
            cmd.png("96sORFstructurecomparison/" + structinfo["Filename"] + ".png", dpi=300)
            cmd.delete("struct1")
        else:
            entry += "NA,NA,NA,NA,"
        if os.path.isfile("PrabakaransORFs2/" + seqdirs + "/couplings/" + prefix + "_CouplingScores.csv"):
            couplingfile = open("PrabakaransORFs2/" + seqdirs + "/couplings/" + prefix + "_CouplingScores.csv","r")
            data = couplingfile.readline()
            data = couplingfile.read()
            couplingfile.close()
            data = data.split("\n")
            data = data[:len(data)-1]
            struct1couplings = {}
            struct1couplingsrank = []
            for line in data:
                temp = line.split(",")
                if(temp[0] < temp[2]):
                    struct1couplings[str(temp[0]) + temp[1] + str(temp[2]) + temp[3]] = temp[5]
                    struct1couplingsrank.append(str(temp[0]) + temp[1] + str(temp[2]) + temp[3])
                else:
                    struct1couplings[str(temp[2]) + temp[3] + str(temp[0]) + temp[1]] = temp[5]
                    struct1couplingsrank.append(str(temp[2]) + temp[3] + str(temp[0]) + temp[1])
        elif os.path.isfile("PrabakaransORFs2/" + seqdirs + "ev_couplings/" + prefix + "_CouplingScores.csv"):
            couplingfile = open("PrabakaransORFs2/" + seqdirs + "/ev_couplings/" + prefix + "_CouplingScores.csv","r")
            data = couplingfile.read()
            couplingfile.close()
            data = data.split("\n")
            data = data[:len(data)-1]
            struct1couplings = {}
            struct1couplingsrank = []
            for line in data:
                temp = line.split(",")
                if(temp[0] < temp[1]):
                    struct1couplings[str(temp[0]) + temp[10] + str(temp[1]) + temp[11]] = temp[2]
                    struct1couplingsrank.append(str(temp[0]) + temp[10] + str(temp[1]) + temp[11])
                else:
                    struct1couplings[str(temp[1]) + temp[11] + str(temp[0]) + temp[10]] = temp[2]
                    struct1couplingsrank.append(str(temp[1]) + temp[11] + str(temp[0]) + temp[10])

        if os.path.isfile("PrabakaransORFs2refseq" + "/" + seqdirs + "/ncORF" + seqdirs[4:] + ".faa"):
            prefix = "ncORF" + seqdirs[4:]
        else:
            prefix = seqdirs
        if os.path.isfile("PrabakaransORFs2refseq" + "/" + seqdirs + "/fold/" + prefix + "_ranking.csv"):
            structinfo = getTopStructure("PrabakaransORFs2refseq" + "/" + seqdirs + "/fold/" + prefix + "_ranking.csv")
            entry += structinfo["Filename"] + "," + structinfo["Score"] + ",refseq6transcriptframe,local,"
            struct2 = "PrabakaransORFs2refseq" + "/" + seqdirs + "/fold/" + structinfo["Filename"]
            copyfile(struct2,"/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2refseqTopPDBEVfoldranking/" + structinfo["Filename"])
            cmd.load(struct2,"struct2")
            cmd.hide("lines","all")
            cmd.show("cartoon","all")
            cmd.color("cyan","struct2")
            cmd.orient("struct2")
            cmd.png("96sORFstructurecomparison/" + structinfo["Filename"] + ".png", dpi=300)
            cmd.delete("struct2")
        elif os.path.isfile("PrabakaransORFs2refseq" + "/" + seqdirs + "/structure_outputs/" + prefix + "_alphabeta_ranking.txt"):
            structinfo = getTopStructure("PrabakaransORFs2refseq" + "/" + seqdirs + "/structure_outputs/" + prefix + "_alphabeta_ranking.txt")
            entry += structinfo["Filename"] + "," + structinfo["Score"] + ",refseq6transcriptframe,server,"
            struct2 = "PrabakaransORFs2refseq" + "/" + seqdirs + "/structure_outputs/" + structinfo["Filename"] + ".pdb"
            copyfile(struct2,"/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2refseqTopPDBEVfoldranking/" + structinfo["Filename"] + ".pdb")
            cmd.load(struct2,"struct2")
            cmd.hide("lines","all")
            cmd.show("cartoon","all")
            cmd.color("cyan","struct2")
            cmd.orient("struct2")
            cmd.png("96sORFstructurecomparison/" + structinfo["Filename"] + ".png", dpi=300)
            cmd.delete("struct2")
        else:
            entry += "NA,NA,NA,NA,"
        if os.path.isfile("PrabakaransORFs2refseq/" + seqdirs + "/couplings/" + prefix + "_CouplingScores.csv"):
            couplingfile = open("PrabakaransORFs2refseq/" + seqdirs + "/couplings/" + prefix + "_CouplingScores.csv","r")
            data = couplingfile.readline()
            data = couplingfile.read()
            couplingfile.close()
            data = data.split("\n")
            data = data[:len(data)-1]
            struct2couplings = {}
            struct2couplingsrank = []
            for line in data:
                temp = line.split(",")
                if(temp[0] < temp[2]):
                    struct2couplings[str(temp[0]) + temp[1] + str(temp[2]) + temp[3]] = temp[5]
                    struct2couplingsrank.append(str(temp[0]) + temp[1] + str(temp[2]) + temp[3])
                else:
                    struct2couplings[str(temp[2]) + temp[3] + str(temp[0]) + temp[1]] = temp[5]
                    struct2couplingsrank.append(str(temp[2]) + temp[3] + str(temp[0]) + temp[1])
        elif os.path.isfile("PrabakaransORFs2refseq/" + seqdirs + "ev_couplings/" + prefix + "_CouplingScores.csv"):
            couplingfile = open("PrabakaransORFs2refseq/" + seqdirs + "/ev_couplings/" + prefix + "_CouplingScores.csv","r")
            data = couplingfile.read()
            couplingfile.close()
            data = data.split("\n")
            data = data[:len(data)-1]
            struct2couplings = {}
            struct2couplingsrank = []
            for line in data:
                temp = line.split(",")
                if(temp[0] < temp[1]):
                    struct1couplings[str(temp[0]) + temp[10] + str(temp[1]) + temp[11]] = temp[2]
                    struct1couplingsrank.append(str(temp[0]) + temp[10] + str(temp[1]) + temp[11])
                else:
                    struct1couplings[str(temp[1]) + temp[11] + str(temp[0]) + temp[10]] = temp[2]
                    struct1couplingsrank.append(str(temp[1]) + temp[11] + str(temp[0]) + temp[10])

        if struct1 != None and struct2 != None:
            struct1pymol = cmd.load(struct1, "struct1")
            struct2pymol = cmd.load(struct2, "struct2")
            cmd.hide("lines","all")
            cmd.show("cartoon","all")
            cmd.color("green","struct1")
            cmd.color("cyan","struct2")
            cmd.orient("all")
            sup = cmd.align("struct1","struct2")
            cmd.png("96sORFstructurecomparison/" + seqdirs + ".png", dpi=300)
            print(sup[0])
            entry += str(sup[0]) + ","
            cmd.delete("struct1")
            cmd.delete("struct2")
        else:
            entry += "NA,"
        if struct1couplings and struct2couplings != None:
            couplingsdiff = 0
            couplingsrankdiff = 0
            count = 0
            skipped = 0
            for key, value in struct1couplings.items():
                if key in struct2couplings:
                    couplingsdiff += abs(float(value) - float(struct2couplings[key]))
                    couplingsrankdiff += abs(struct1couplingsrank.index(key) - struct2couplingsrank.index(key))
                    count += 1
                else:
                    skipped += 1
            couplingsdiff /= count
            couplingsrankdiff /= count
            entry += str(couplingsdiff) + "," + str(couplingsrankdiff) + ","
            couplingsrankdiff = float(couplingsrankdiff) / count
            print(couplingsrankdiff)
            print(count)
            entry += str(couplingsrankdiff) + "," + str(skipped) + "\n"
        else:
            entry += "NA,NA,NA,NA\n"
        f.write(entry)
f.close()

