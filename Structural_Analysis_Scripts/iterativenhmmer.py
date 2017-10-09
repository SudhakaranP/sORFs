#!/usr/bin/python
import sys
import os
import subprocess

PATH_TO_DATA = "/media/david/genomicsdrive/NovelPeptideCharacterisation/data"
#PATH_TO_OUTPUT = "/media/david/genomicsdrive/NovelPeptideCharacterisation/HS_genbank_GRCh38_v1"
PATH_TO_OUTPUT = "/media/david/genomicsdrive/NovelPeptideCharacterisation/PrabakaransORFs"

def backtoprotein(seqdir):
    args = ["transeq", "-sequence", seqdir + "-nhmmer-final.a2m",
            "-outseq", seqdir + "-nhmmer-final.faa"]
    p = subprocess.check_output(args)
    args = ["phmmer", "-A", seqdir + "-phmmer-final.sto",
            "-o", seqdir + "-phmmer-out.txt",
            seqdir + ".faa", seqdir + "-nhmmer-final.faa"]
    p = subprocess.check_output(args)
    args = ["sreformat", "a2m", seqdir + "-phmmer-final.sto"]
    p = subprocess.checkoutput(args)
    f = open(seqdir + "-phmmer-final.a2m", "w")
    f.write(p)
    f.close()

def CheckIteration(seqdir, i, seqs, final):
    seqadded = 0
    seqdropped = 0
    #Convert alignemnt to fasta format
    args = ["sreformat", "fasta", seqdir + "-nhmmer-" + str(i-1) + ".sto"]
    p = subprocess.check_output(args)
    tmp = open(seqdir + "-tmp.fna","w")
    tmp.write(p)
    tmp.close()
    #Translate to protein
    args = ["transeq","-sequence",seqdir + "-tmp.fna", "-outseq", seqdir + "-tmp.faa"]
    p = subprocess.check_output(args)
    #Check for significant alignment
    args = ["phmmer","-A",seqdir + "-phmmer-" + str(i-1) + ".sto",
            seqdir + ".faa", seqdir + "-tmp.faa"]
    p = subprocess.check_output(args)
    if os.stat(seqdir + "-phmmer-" + str(i-1) + ".sto").st_size == 0:
        tmp = open(seqdir + "-nhmmer-final.a2m","w")
        tmp.write("")
        tmp.close()
        print("No significant hits found")
        return
    else:
        args = ["sreformat","fasta",seqdir + "-phmmer-" + str(i-1) + ".sto"]
        p = subprocess.check_output(args)
        tmp = open(seqdir + "-phmmer-" + str(i-1) + ".faa", "w")
        tmp.write(p)
        tmp.close()
    #Create list of hits in nucleotide form
    args = ["backtranseq","-sequence",seqdir + "-phmmer-" + str(i-1) + ".faa",
            "-outfile",seqdir + "-phmmer-" + str(i-1) + ".fna"]
    p = subprocess.check_output(args)
    tmp = open(seqdir + ".fna","r")
    originalseq = tmp.read()
    tmp.close()
    tmp = open(seqdir + "-phmmer-" + str(i-1) + ".fna","r")
    hits = tmp.read()
    tmp.close()
    newseqs = originalseq + hits
    tmp = open(seqdir + "-phmmer-" + str(i-1) + ".fna","w")
    tmp.write(newseqs)
    tmp.close()
    #Check for convergence
    lines = newseqs.split("\n")
    for l in lines:
        if len(l) > 0:
            if l[0] == ">":
                if l not in seqs:
                    seqs.append(l)
                    seqadded += 1
    print("Sequences added: " + str(seqadded))
    for s in seqs:
        if s not in lines:
            seqs.remove(s)
            seqdropped += 1
    print("Sequences dropped: " + str(seqadded))
    #Create alignment file for hmmbuild
    args = ["nhmmer", "-A", seqdir + "-nhmmer-" + str(i-1) + ".sto",
            seqdir + ".fna", seqdir + "-phmmer-" + str(i-1) + ".fna"]
    p = subprocess.check_output(args)
    if (seqadded == 0 and seqdropped == 0) or final:
        args = ["sreformat", "a2m", seqdir + "-nhmmer-" + str(i-1) + ".sto"]
        p = subprocess.check_output(args)
        tmp = open(seqdir + "-nhmmer-final.a2m","w")
        tmp.write(p)
        tmp.close()
        if not final:
            print("Search has converged")
        return False
    return True

def itnhmmer(path, seqdir, it, override):
    os.chdir(path)
    seqs = []
    startpoint = 0
    if not override:
        skipiter = False
        while startpoint < it and os.path.isfile(seqdir + "-nhmmer-" + str(startpoint) + ".hmm"):
            startpoint += 1
            skipiter = True
        if skipiter:
            startpoint -= 1
    i = 0
    for i in range(startpoint, it):
        print "Round " + str(i)
        if i == 0:
            args = ["hmmbuild", seqdir + "-nhmmer-" + str(i) + ".hmm",
                    seqdir + ".fna"]
            tmp = open(seqdir + ".fna","r")
            seqs.append(tmp.readline())
            tmp.close()
            runnhmmer = True
            #Build hmm file
            p = subprocess.check_output(args)
            print p
        else:
            runnhmmer = CheckIteration(seqdir,i,seqs, False)
            #Build hmm file
            args = ["hmmbuild", seqdir + "-nhmmer-" + str(i) + ".hmm",
                seqdir + "-nhmmer-" + str(i-1) + ".sto"]
            p = subprocess.check_output(args)
            print p
        if runnhmmer:
            #Execute nhmmer alignment
            args = ["nhmmer",
                "-o", path + "/" + seqdir + "-nhmmer-out-" + str(i) + ".txt",
                "-A", path + "/" + seqdir + "-nhmmer-" + str(i) + ".sto",
                seqdir + "-nhmmer-" + str(i) + ".hmm",
                PATH_TO_DATA + "/nt"]
            p = subprocess.check_output(args)
        else:
            break
    CheckIteration(seqdir,i,seqs, True)
    backtoprotein(seqdir)

def executenhmmer(path, seqdir, override):
    print "Processing " + seqdir
    if not override:
        if not os.path.isfile(path + "/" + seqdir + "-nhmmer-final.a2m"):
            itnhmmer(path, seqdir, 5, override)
    else:
        itnhmmer(path, seqdir, 5, override)

argpos = 1
override = False

if len(sys.argv) == 1:
    print "Usage:"
    print "bin/iterativenhmmer.py <-o> <-all/-list/-lim> <val>"
    print "-o \t forces nhmmer to recompute even if files exists"
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
            executenhmmer(root + "/" + seqdir, seqdir, override)
elif sys.argv[argpos] == '-lim':
    counter = 0
    lim = int(sys.argv[argpos + 1])
    for root, dirs, files in os.walk(PATH_TO_OUTPUT):
        for seqdir in dirs:
            executenhmmer(root + "/" + seqdir, seqdir, override)
            counter = counter + 1
            if counter == lim:
                        sys.exit()
elif sys.argv[argpos] == '-list':
    dirs = sys.argv[argpos + 1:]
    for seqdir in dirs:
        executenhmmer(PATH_TO_OUTPUT + "/" + seqdir, seqdir, override)
else:
    print "Usage"
    print "bin/iterativenhmmer.py <-o> <-all/-list/-lim> <val>"
    print "-o \t forces jackhmmer to recompute even if files exists"
    print "-all \t loop over all directories in output file"
    print "-list \t loops over val which must be a list of directory names separated by spaces"
    print "-lim \t loops over a certain number of directories provided as val before terminating"
