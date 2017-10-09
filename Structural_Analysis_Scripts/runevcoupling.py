#!usr/bin/python
import sys
from evcouplings.utils import read_config_file
from evcouplings.utils.pipeline import execute

CONFIGPATH = "EVcouplings-master/config/protein-struct-pred-config.txt"

if len(sys.argv) == 1:
    print("runevcoupling.py <outputdir> <sequenceid>")
else:
    outputdir = sys.argv[1]
    seqid = sys.argv[2]
    config = read_config_file(CONFIGPATH, preserve_order = True)
    config["global"]["prefix"] = outputdir
    config["global"]["sequence_id"] = seqid
    config["global"]["sequence_file"] = outputdir + "/" + seqid + ".faa"
    outcfg = execute(**config)
