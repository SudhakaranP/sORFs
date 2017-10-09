#!/usr/bin/python
import os
import sys
from pymol import *

dirpath = "/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs2refseqTopPDBEVfoldranking/"

f = open("PrabakaransORFs2refseqTopPDBEVfoldrmsdmat.csv","w")
f.write(",".join(os.listdir(dirpath)) + "\n")
for pdbfile in os.listdir(dirpath):
    entry = ""
    cmd.load(dirpath + pdbfile, "curpdb")
    for secondfile in os.listdir(dirpath):
        if secondfile == pdbfile:
            entry += "0,"
        else:
            cmd.load(dirpath + secondfile,"testpdb")
            sup = cmd.super("curpdb","testpdb")
            print(sup[0])
            entry += str(sup[0]) + ","
            cmd.delete("testpdb")
    entry = entry[:len(entry)-1] + "\n"
    cmd.delete("curpdb")
    f.write(entry)
f.close()