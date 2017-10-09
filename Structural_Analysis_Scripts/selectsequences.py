#!/usr/bin/python
import os
import sys
import numpy as np
import datetime as dt

PATH_TO_DATA = "/media/david/Seagate Backup Plus Drive/Summer 2017 Genetics Project/data"
PATH_TO_OUTPUT = "/media/david/Seagate Backup Plus Drive/Summer 2017 Genetics Project/output"
q = 21
l = 0.5
x = 0.7
iterthreshold = 8


def getseqaslist(line):
    temp = line.split(" ")
    seqid = temp[0]
    seq = temp[len(temp) - 1]
    seq = seq[:len(seq) - 1]
    temp = []
    for a in seq:
        temp.append(a)
    return (seqid, temp)

def countdiscount(line, indices, char):
    count = 0
    for i in range(0, len(line)):
        if i not in indices:
            if line[i] == char:
                count = count + 1
    return count


def loadfileintomat(line):
    overallalignment = 0
    with open(line) as alignmentfile:
        start = False
        msamatlist = []
        seqidlist = []
        msamat = None
        seqids = []
        seqscores = []
        overallalignment = []
        seqlen = 0
        for line in alignmentfile:
            if not start:
                if not (line[0] == "\n" or line[0] == "#" or line[0] == "/"):
                    start = True
                    parsedline = getseqaslist(line)
                    msamat = np.array(parsedline[1])
                    seqids.append(parsedline[0])
                    seqlen = msamat.shape[0]
            else:
                if line[0] == "/" or line[0] == "\n":
                    pass
                elif line[0] != "#":
                    parsedline = getseqaslist(line)
                    newseq = np.array(parsedline[1])
                    msamat = np.vstack((msamat, newseq))
                    seqids.append(parsedline[0])
                elif line[0:4] == "#=GR":
                    #if alignment is less than 50% of query discard
                    parsedline = getseqaslist(line)
                    seqscores.append(parsedline[1])
                elif line[0:7] == "#=GC RF":
                    #remove columns with excessive gaps
                    print(msamat.shape)
                    colstoremove = []
                    if len(msamat.shape) > 1:
                        for i in range(0, msamat.shape[1]):
                            gapcount = 0
                            for row in msamat:
                                if row[i] == "-":
                                    gapcount = gapcount + 1
                                else:
                                    gapcount = gapcount - 1
                            if gapcount > 0:
                                colstoremove.append(i)
                        msamat = np.delete(msamat, colstoremove, 1)
                    print("Removed cols:")
                    print msamat.shape
                    #return consensus percentage for checking overall alignment
                    temp = getseqaslist(line)
                    seqlen = seqlen - len(colstoremove)
                    alignmentcount = countdiscount(temp[1],colstoremove,"x")
                    overallalignment.append(float(alignmentcount) / seqlen)
                    #Remove sequences with less than 50% alignment to reference
                    rowstoremove = []
                    for i in range(0, len(seqscores)):
                        scorecount = countdiscount(seqscores[i],colstoremove, ".")
                        if scorecount > 0.5 * seqlen:
                            rowstoremove.append(i)
                    msamat = np.delete(msamat, rowstoremove, 0)
                    print("Removed rows:")
                    print(msamat.shape)
                    msamatlist.append(msamat)
                    adjustedseqids = []
                    for i in range(0, len(seqidlist)):
                        if i not in rowstoremove:
                            adjustedseqids.append(seqids[i])
                    seqidlist.append(seqids)
                    msamat = None
                    start = False
        return (msamatlist, overallalignment, seqidlist)


def calculatekm(msamat):
    dim = msamat.shape
    if len(dim) == 1:
        nrow = 1
        ncol = dim[0]
    else:
        nrow = dim[0]
        ncol = dim[1]
    kmvec = []
    for n in range(0, nrow):
        neighbourcount = 0
        for m in range(0, nrow):
            tempsum = 0
            for i in range(0, ncol):
                if nrow > 1:
                    if msamat[n, i] == msamat[m, i]:
                        tempsum = tempsum + 1
                else:
                    tempsum = ncol
                    break
            tempsum = tempsum - (x * ncol)
            if tempsum >= 0:
                neighbourcount = neighbourcount + 1
        kmvec.append(float(neighbourcount))
    return np.array(kmvec)

for seqdir in os.listdir(PATH_TO_OUTPUT):
    ts = dt.datetime.now()
    output = open(PATH_TO_OUTPUT + "/selectedseq-x" + str(x * 100) + "-iter" +
            str(iterthreshold) + "-" + ts.strftime("%Y-%m-%d") + ".txt", "w")
    if os.path.isfile(PATH_TO_OUTPUT + "/" + seqdir + "/" +
                        seqdir + "-" + str(iterthreshold) + ".sto"):
        print "Processing " + seqdir
        data = loadfileintomat(PATH_TO_OUTPUT + "/" + seqdir + "/"
                + seqdir + "-" + str(iterthreshold) + ".sto")
        msamatlist = data[0]
        overallalignmentlist = data[1]
        for i in range(0, len(msamatlist)):
            if overallalignmentlist[i] > 0.95:
                dim = msamatlist[i].shape
                if len(dim) == 1:
                    nrow = 1
                    ncol = dim[0]
                else:
                    nrow = dim[0]
                    ncol = dim[1]
                kmvec = calculatekm(msamatlist[i])
                print(kmvec)
                mef = sum(kmvec ** - 1)
                print("Mef per residue: " + str(mef) + " / " + str(ncol) + "=" + str(float(mef / ncol)))
                print("Overall alignment: " + str(overallalignmentlist[i]))
                if float(mef / ncol) > 5:
                    a2mfile = open(PATH_TO_OUTPUT + "/" + seqdir + "/" + seqdir +
                            "_alignment.a2m","w")
                    for j in range(0, ncol):
                        a2mfile.write(">" + data[2][i][j] + "\n")
                        a2mfile.write(msamatlist[i][j, ].tostring() + "\n")
                    a2mfile.close()
                    output.write(seqdir + "\n")
print "Done"
output.close()
sys.exit()



