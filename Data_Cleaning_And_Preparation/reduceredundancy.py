#!/usr/bin/python
import sys

PATH_TO_DATA = "/media/david/Seagate Backup Plus Drive/Summer 2017 Genetics Project/data"

def removeduplicates(seqs,seqids):
    count += 1
    print("Stack Size: " + str(count))
    if seqs == []:
        return ([],[])
    elif len(seqs) == 1:
        return (seqs,seqids)
    else:
        pivot = int(len(seqs) / 2)
        return mergelists(removeduplicates(seqs[:pivot],seqids[:pivot]),
                    removeduplicates(seqs[pivot:],seqids[pivot:]),count)
def mergelists(t1,t2,count):
    print("Resolving stack: " + str(count))
    for i in range(0,len(t2[0])):
        if t2[0][i] not in t1[0]:
            t1[0].append(t2[0][i])
            t1[1].append(t2[1][i])
    return t1

if len(sys.argv) == 1:
    print("Please include a fasta file")
else:
    seqids = [] #Hold sequence IDs
    seqs = [] #Hold sequence data
    curseqid = ""
    curseq = ""
    redundantcount = 0 #Count total sequences in original file
    count = 0 #Count final sequences after reduction
    seqcount = 0 #Count number of sequences in seqs
    partcount = 0 #Count current part number of database split
    with open(PATH_TO_DATA + "/" + sys.argv[1], "r") as seqfile:
        header = ""
        for line in seqfile:
            if line[0] == ">":
                if curseqid != "":
                    print(seqcount)
                    seqs.append(curseq)
                    seqids.append(curseqid)
                    curseq = ""
                    seqcount += 1
                curseqid = line
                redundantcount += 1
            else:
                curseq += line

            if seqcount == 1000000:
                seqcount = 0
                partcount += 1
                print("Processing part" + str(partcount))
                results = removeduplicates(seqs,seqids)
                f = open(PATH_TO_DATA + "/part-" + str(partcount) + "-" +
                            sys.argv[1], "w")
                for i in range(0,len(results[0])):
                    f.write(results[1][i] + results[0][i])
                f.close()
                seqs = []
                seqids = []

    seqs.append(curseq)
    seqids.append(curseqid)
    seqcount = 0
    partcount += 1
    print("Processing part" + str(partcount))
    results = removeduplicates(seqs,seqids)
    f = open(PATH_TO_DATA + "/part-" + str(partcount) + "-" +
                sys.argv[1], "w")
    for i in range(0,len(results[0])):
        f.write(results[1][i] + results[0][i])
    f.close()
    #Begin merging parts

    #nonredundantseqs = removeduplicates(seqs,seqids)
    #print "Writing sequences"
    #f = open(PATH_TO_DATA + "/non-redundant-" + sys.argv[1], "w")
    #for i in range(0,len(nonredundantseqs[0])):
    #    f.write(nonredundantseqs[1][i] + nonredundantseqs[0][i])
    #f.close()
    #print(str(redundantcount) + " Sequences processed")
    #print(str(count) + " Sequences remain")