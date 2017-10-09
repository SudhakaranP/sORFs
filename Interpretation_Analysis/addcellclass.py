#!/usr/bin/python

f = open("3604_btclassification.txt","r")
btclass = f.read()
f.close()
f = open("data/peptides_categories_table.tsv","r")
neurongeneexp = f.read()
f.close()
f = open("exonicgenelistncbiid.csv","r")
exonicgenelist = f.read()
f.close()
f = open("data/PrabakaransORFsTranslationEvidence.faa","r")
fasta = f.read()
f.close()

btclass = btclass.split("\n")
neurongeneexp = neurongeneexp.split("\n")
exonicgenelist = exonicgenelist.split("\n")
fasta = fasta.split("\n")


neurongenelist = []
for candidate in neurongeneexp:
    if len(candidate) > 0:
        temp = candidate.split("\t")
        if temp[14].strip() != "-":
            for gene in exonicgenelist:
                if len(gene) > 0:
                    temp2 = gene.split(",")
                    if temp2[1].upper() == temp[14].strip().upper():
                        neurongenelist.append(temp2)
T = 0
B = 0
BN = 0
TN = 0
BT = 0
BTN = 0

f = open("data/PrabakaransORFsTranslationEvidence.faa","w")
for line in fasta:
    if len(line) > 0:
        if line[0] == ">":
            temp = line.split("|")
            term = ""
            for sorf in btclass:
                if len(sorf) > 0:
                    sorfterms = sorf.split("\t")
                    if temp[1] == sorfterms[0]:
                        if sorfterms[1] == "BCELL":
                            term += "B"
                        if sorfterms[2] == "TCELL":
                            term += "T"
                        break
            for gene in neurongenelist:
                if temp[4] == gene[2]:
                    if (temp[5] < gene[3] and temp[6] > gene[3]) or (temp[5] > gene[3] and temp[5] < gene[4]):
                        term += "N"
                        break
            if temp[8] == "exonic":
                if term == "B":
                    B += 1
                elif term == "T":
                    T += 1
                elif term == "BT":
                    BT += 1
                elif term == "BN":
                    BN += 1
                elif term == "TN":
                    TN += 1
                elif term == "BTN":
                    BTN += 1

            if len(temp) == 9:
                line += "|" + term
            else:
                temp2 = line.split("|")
                temp2[9] == term
                line = "|".join(temp2)
    print(line)
    f.write(line + "\n")
f.close()

print("B:" + str(B))
print("T:" + str(T))
print("BN:" + str(BN))
print("TN:" + str(TN))
print("BT:" + str(BT))
print("BTN:" + str(BTN))
print("Total:" + str(B + T + BT + TN + BN +BTN))
print(len(neurongenelist))