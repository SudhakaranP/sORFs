#!/usr/bin/python
import sys
import os
import subprocess

PATH_TO_DATA = "/media/david/genomicsdrive/NovelPeptideCharacterisation/data"
PATH_TO_OUTPUT = "/media/david/genomicsdrive/NovelPeptideCharacterisation/HS_genbank_hg19_v1"


def executejackhmmer(path, override):
    args = ["jackhmmer", "--cpu", "4",
            "--chkhmm", path + "/" + seqdir + "uniref",
            "--chkali", path + "/" + seqdir + "uniref",
            "--tblout", path + "/" + seqdir + "_jackhmmer_uniref.txt",
            "-A", path + "/" + seqdir + "_alignment_uniref.sto",
            "-o", path + "/" + seqdir + "_jackhmmer_verbose_uniref.txt",
            "-N", "11",
            path + "/" + seqdir + ".fa",
            PATH_TO_DATA + "/" + "uniref100.fasta"]
    if not override:
        if not os.path.isfile(path + "/" + seqdir + "_jackhmmer_uniref.txt"):
            p = subprocess.check_output(args)
            print p
    else:
        p = subprocess.check_output(args)
        print p

argpos = 1
override = False

if len(sys.argv) == 1:
    print "Usage:"
    print "bin/jackhmmer.py <-o> <-all/-list/-lim> <val>"
    print "-o \t forces jackhmmer to recompute even if files exists"
    print "-all \t loop over all directories in output file"
    print "-list \t loops over val which must be a list of directory names separated by spaces"
    print "-lim \t loops over a certain number of directories provided as val before terminating"
    sys.exit()
elif sys.argv[argpos] == "-o":
    argpos = argpos + 1
    override = True

if sys.argv[argpos] == '-all':
    for root, dirs, files in os.walk(PATH_TO_OUTPUT):
        for seqdir in dirs:
            print("Processing " + seqdir)
            executejackhmmer(root + "/" + seqdir, override)
elif sys.argv[argpos] == '-lim':
    counter = 0
    lim = int(sys.argv[argpos + 1])
    for root, dirs, files in os.walk(PATH_TO_OUTPUT):
        for seqdir in dirs:
            print("Processing " + seqdir)
            executejackhmmer(root + "/" + seqdir, override)
            counter = counter + 1
            if counter == lim:
                        sys.exit()
elif sys.argv[argpos] == '-list':
    dirs = sys.argv[argpos + 1:]
    for seqdir in dirs:
        print("Processing " + seqdir)
        executejackhmmer(PATH_TO_OUTPUT + "/" + seqdir, override)
else:
    print "Usage"
    print "bin/jackhmmer.py <-o> <-all/-list/-lim> <val>"
    print "-o \t forces jackhmmer to recompute even if files exists"
    print "-all \t loop over all directories in output file"
    print "-list \t loops over val which must be a list of directory names separated by spaces"
    print "-lim \t loops over a certain number of directories provided as val before terminating"

