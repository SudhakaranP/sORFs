#!usr/bin/python
import sys
import random
import math

def perturb(pose,scorefxn):
    res = random.randint(1,len(pose.sequence()))
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
                pose.set_phi(res,pose.phi(res - delta))
            else:
                pose.set_psi(res,pose.phi(res - delta))
            return pose

def nextres(pose, res, lowpose):
    lowestenergy = -1
    newres = pose_from_sequence(res)
    pose = insert_pose_into_pose(pose, newres, len(pose.sequence()))
    for range in (0,500000):
        pose = perturb(pose,scorefxn)
        if lowestenergy < 0:
            lowestenergy = scorefxn(pose)
        else:
            curenergy = scorefxn(pose)
            if lowestenergy > curenergy:
                lowpose.assign(pose)
                lowestenergy = curenergy
    print "Lowest Energy: " + str(lowestenergy)
    print "Current Energy: " + str(scorefxn(pose))
    lowpose.dump_pdb(lowpose.sequence() + ".pdb")

if len(sys.argv) == 1:
    print "montecarlofold.py <amino acid sequence>"
else:
    from pyrosetta import *
    from rosetta.core.scoring import *
    from rosetta.protocols.grafting import *
    init()
    if len(sys.argv[1]) > 2:
        respos = 1
        pose = pose_from_sequence(sys.argv[1][0:respos])
        lowpose = pose_from_sequence('')
        scorefxn = ScoreFunction()
        scorefxn.set_weight(fa_atr, 1.0)
        scorefxn.set_weight(fa_rep, 1.0)
        scorefxn.set_weight(hbond_sr_bb, 1.0)
        scorefxn.set_weight(hbond_lr_bb, 1.0)
        scorefxn.set_weight(hbond_bb_sc, 1.0)
        scorefxn.set_weight(hbond_sc, 1.0)
        while respos < len(sys.argv[1]):
            print pose.sequence()
            nextres(pose, sys.argv[1][respos], lowpose)
            respos += 1