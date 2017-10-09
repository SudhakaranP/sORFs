#!/usr/bin/python
import os

f = open("/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2-PrabakaransORFs2refseq-compare.csv","r")
header = f.readline()
data = f.read()
f.close()
data = data.split("\n")
data = data[:len(data)-1]
header = header.split(",")
header.remove("database1")
header.remove("foldingsoftware1")
header.remove("database2")
header.remove("foldingsoftware2")
header.insert(5, "superimpose")

f = open("/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2-PrabakaransORFs2refseq-compare.html","w")
f.write("<html><body><table border=1 border-collapse>")
f.write("<tr><th>Info</th><th>Struct1</th><th>Struct2</th><th>Superimpose</th></tr>")
for line in data:
    fields = line.split(",")
    f.write("<tr>")
    f.write("<td>ID: " + fields[0] + "<br />")
    if fields[10] != "NA":
        f.write("Coupling Prob Dif:" + str(round(float(fields[10]),3)) + "<br />")
    else:
        f.write("Coupling Prob Dif: NA<br />")
    if fields[11] != "NA":
        f.write("Coupling Rank Dif: " + str(round(float(fields[11]),3)) + "<br />")
    else:
        f.write("Coupling Rank Dif: NA<br />")
    if fields[12] != "NA":
        f.write("Normalised Coupling Rank Dif: " + str(round(float(fields[12]),3)) + "<br />")
    else:
        f.write("Coupling Rank Dif: NA<br />")
    if fields[13] != "NA":
        f.write("Skipped Couplings: " + str(round(float(fields[13]),3)) + "<br />")
    else:
        f.write("Skipped Couplings: NA</td>")

    if fields[1] != "NA":
        f.write("<td><img src='/media/david/genomicsdrive/NovelPeptideCharacterisation/96sORFstructurecomparison/" + fields[1] + ".png' width=350><br />")
        if fields[2] != "NA":
            f.write("alphabeta score: " + str(round(float(fields[2]),3)) + "</td>")
        else:
            f.write("alphabeta score: NA</td>")
    else:
        f.write("<td>NA</td>")

    if fields[5] != "NA":
        f.write("<td><img src='/media/david/genomicsdrive/NovelPeptideCharacterisation/96sORFstructurecomparison/" + fields[5] + ".png' width=350><br />")
        if fields[6] != "NA":
            f.write("alphabeta score: " + str(round(float(fields[6]),3)) + "</td>")
        else:
            f.write("alphabeta score: NA</td>")
    else:
        f.write("<td>NA</td>")

    if fields[1] != "NA" and fields[5] != "NA":
        f.write("<td><img src='/media/david/genomicsdrive/NovelPeptideCharacterisation/96sORFstructurecomparison/" + fields[0] + ".png' width=350><br />")
        if fields[9] != "NA":
            f.write("RMSD: " + str(round(float(fields[9]),3)) + "</td>")
        else:
            f.write("RMSD: NA</td>")
    else:
        f.write("<td>NA</td>")

    f.write("</tr>")
f.write("</table></html>")