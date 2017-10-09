#!usr/bin/python
import sys
import random
import math

def perturb(pose,scorefxn,validpositions):
    res = validpositions[random.randint(0,len(validpositions)-1)]
    preve = scorefxn(pose)
    delta = random.gauss(0,25)
    phiorsi = random.randint(0,1)
    if phiorsi == 1:
        pose.set_phi(res,pose.phi(res) + delta)
    else:
        pose.set_psi(res,pose.psi(res) + delta)

    poste = scorefxn(pose)
    deltae = poste - preve
    if deltae < 0:
        return pose
    else:
        p = math.exp(-deltae)
        if p > 0.5:
            return pose
        else:
            if phiorsi == 1:
                pose.set_phi(res,pose.phi(res)- delta)
            else:
                pose.set_psi(res,pose.phi(res) - delta)
            return pose

if len(sys.argv) == 1:
    print "montecarlofold.py <seq file> <sec struct>"
else:
    from pyrosetta import *
    from rosetta.core.scoring import *
    init()
    f = open(sys.argv[1],"r")
    seq = ""
    for line in f:
        if line[0] != ">":
            seq += line.strip()
    f.close()
    f = open(sys.argv[2], "r")
    secstruct = ""
    for line in f:
        if len(line) > 5:
            if line[0:4] == "Pred":
                temp = line.split(" ")
                secstruct += temp[1].strip()
    f.close()
    pose = pose_from_sequence(seq);
    lowpose = pose_from_sequence("")
    scorefxn = ScoreFunction()
    scorefxn.set_weight(fa_atr, 1.0)
    scorefxn.set_weight(fa_rep, 1.0)
    scorefxn.set_weight(hbond_sr_bb, 1.0)
    scorefxn.set_weight(hbond_lr_bb, 1.0)
    scorefxn.set_weight(hbond_bb_sc, 1.0)
    scorefxn.set_weight(hbond_sc, 1.0)
    lowestenergy = -1

    validpositions = []
    print scorefxn(pose)
    for i in range(0,len(secstruct)-1):
        if secstruct[i] == "H":
            pose.set_phi(i,-64)
            pose.set_psi(i,-41)
        elif secstruct[i] == "E":
            pose.set_phi(i,-135)
            pose.set_psi(i,135)
        else:
            validpositions.append(i+1)
    print scorefxn(pose)

    for i in range(0,5000 * len(pose.sequence())):
        pose = perturb(pose,scorefxn,validpositions)
        if lowestenergy < 0:
            lowestenergy = scorefxn(pose)
        else:
            curenergy = scorefxn(pose)
            if lowestenergy > curenergy:
                lowpose.assign(pose)
                lowestenergy = curenergy
    print "Lowest Energy: " + str(lowestenergy)
    print "Current Energy: " + str(scorefxn(pose))
    lowpose.dump_pdb(sys.argv[1].split(".")[0] + ".pdb")