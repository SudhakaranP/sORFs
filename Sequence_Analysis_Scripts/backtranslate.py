#!/usr/bin/python
import os
import sys
import subprocess

PATH_TO_OUTPUT = "/media/david/genomicsdrive/NovelPeptideCharacterisation/"

if len(sys.argv) == 1:
    print("backtranslate.py <Dir>")
else:
    for root, dirs, files in os.walk(PATH_TO_OUTPUT + sys.argv[1]):
        for seqdir in dirs:
            seqfilepath = root + "/" + seqdir + "/" + seqdir + ".faa"
            translatedfilepath = root + "/" + seqdir + "/" + seqdir + ".fna"
            if os.path.isfile(seqfilepath):
                args = ["backtranseq",
                    "-sequence", seqfilepath,
                    "-outfile", translatedfilepath]
                p = subprocess.check_output(args)
                print p
