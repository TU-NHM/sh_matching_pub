#!/usr/bin/python
import argparse
import logging
import os
import re
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to remove sequences with more than X ambiguous bases)")
parser.add_argument("run_id", help="Need run id in numeric format!")
parser.add_argument("allowed_number", help="Need allowed number of ambiguous bases in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
allowed_number = args.allowed_number
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)
if not allowed_number.isdigit():
    raise ValueError("Allowed number of ambiguous bases is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile = user_dir / "seqs_out_chim.fasta"
outfile = user_dir / "iupac_out.fasta"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

infile_hash1 = {}

# open excluded seq file
with open(ex_file, "a") as ex, open(outfile, "w") as o, open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        seq = str(record.seq).upper()
        sequencelength = len(seq)
        regexPattern = re.compile('[YRSWKMBDHVN]')
        listOfmatches = regexPattern.findall(seq)
        numberofIUPACs = len(listOfmatches)
        percent = (numberofIUPACs / sequencelength) * 100
        if numberofIUPACs > int(allowed_number):
            n_search = re.search("^([A-Z]+?)(N{6,})(.*)$", seq)
            if n_search:
                part_1 = n_search.group(1)
                part_2 = n_search.group(3)
                sequencelength2 = len(part_1)
                if sequencelength2 >= 300:
                    listOfmatches2 = regexPattern.findall(part_1)
                    numberofIUPACs2 = len(listOfmatches2)
                    percent2 = (numberofIUPACs2 / sequencelength2) * 100
                    if numberofIUPACs2 > int(allowed_number):
                        logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs}, cut)")
                        ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs}) in cut sequence exceeds the number of allowed ambiguous bases ({allowed_number})\n")
                    else:
                        o.write(f">{name}\n")
                        o.write(f"{part_1}\n")
                else:
                    logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs}, cut but too short to fix)")
                    ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs}) in sequence exceeds the number of allowed ambiguous bases ({allowed_number}). Cut, but too short to fix.\n")
            else:
                logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs})")
                ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs}) in sequence exceeds the number of allowed ambiguous bases ({allowed_number})\n")
        else:
            o.write(f">{name}\n")
            o.write(f"{seq}\n")
