#!/usr/bin/python
import sys

PATH_TO_DATA = "/media/david/Seagate Backup Plus Drive/Summer 2017 Genetics Project/data"

def removeduplicates(seqs,seqids,count):
    count += 1
    print("Stack Size: " + str(count))
    if seqs == []:
        return ([],[])
    elif len(seqs) == 1:
        return (seqs,seqids)
    else:
        pivot = int(len(seqs) / 2)
        return mergelists(removeduplicates(seqs[:pivot],seqids[:pivot],count),
                    removeduplicates(seqs[pivot:],seqids[pivot:],count),count)
def mergelists(t1,t2,count):
    print("Resolving stack: " + str(count))
    for i in range(0,len(t2[0])):
        if t2[0][i] not in t1[0]:
            t1[0].append(t2[0][i])
            t1[1].append(t2[1][i])
        else:
            print("Sequence Removed")
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
                    seqs.append(curseq)
                    seqids.append(curseqid)
                    seqcount += 1
                    curseq = ""
                curseqid = line
                redundantcount += 1
            else:
                curseq += line

    seqs.append(curseq)
    seqids.append(curseqid)
    nonredundantseqs = removeduplicates(seqs,seqids,0)
    print "Writing sequences"
    f = open(PATH_TO_DATA + "/non-redundant-" + sys.argv[1], "w")
    for i in range(0,len(nonredundantseqs[0])):
        f.write(nonredundantseqs[1][i] + nonredundantseqs[0][i])
    f.close()
    print(str(redundantcount) + " Sequences processed")
    print(str(count) + " Sequences remain")