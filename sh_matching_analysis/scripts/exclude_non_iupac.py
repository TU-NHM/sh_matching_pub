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
infile_orig = user_dir / f"source_{run_id}_fastaunique"
outfile = user_dir / "iupac_out_full.fasta"
outfile_vsearch_96 = user_dir / "iupac_out_vsearch_96.fasta"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

# use the original sequence for vsearch 100% clustering and not the one that came out from ITSx extractor
orig_seqs_dict = {}
with open(infile_orig, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        orig_seqs_dict[record.id] = record.seq

infile_hash1 = {}

# open excluded seq file
with open(ex_file, "a") as ex, open(outfile, "w") as o, open(outfile_vsearch_96, "w") as o2, open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        seq = str(record.seq).upper()
        sequencelength = len(seq)
        regexPattern_1 = re.compile('[YRSWKMBDHVN]')
        regexPattern_2 = re.compile('[N]')
        listOfmatches_1 = regexPattern_1.findall(seq)
        listOfmatches_2 = regexPattern_2.findall(seq)
        numberofIUPACs_1 = len(listOfmatches_1)
        numberofIUPACs_2 = len(listOfmatches_2)
        if numberofIUPACs_2 > int(allowed_number) or numberofIUPACs_1 > 16:
            n_search = re.search("^([A-Z]+?)(N{6,})(.*)$", seq)
            if n_search:
                part_1 = n_search.group(1)
                part_2 = n_search.group(3)
                sequencelength2 = len(part_1)
                if sequencelength2 >= 300:
                    listOfmatches2_1 = regexPattern_1.findall(part_1)
                    listOfmatches2_2 = regexPattern_2.findall(part_1)
                    numberofIUPACs2_1 = len(listOfmatches2_1)
                    numberofIUPACs2_2 = len(listOfmatches2_2)
                    if numberofIUPACs2_2 > int(allowed_number) or numberofIUPACs2_1 > 16:
                        logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs_1}, {numberofIUPACs_2}, cut)")
                        ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs_1}, {numberofIUPACs_2}) in cut sequence exceeds the number of allowed ambiguous bases (16, {allowed_number})\n")
                    else:
                        o.write(f">{name}\n")
                        o.write(f"{part_1}\n")
                        o2.write(f">{name}\n")
                        o2.write(f"{orig_seqs_dict[name]}\n")
                else:
                    logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs_1}, {numberofIUPACs_2}, cut but too short to fix)")
                    ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs_1}, {numberofIUPACs_2}) in sequence exceeds the number of allowed ambiguous bases (16, {allowed_number}). Cut, but too short to fix.\n")
            else:
                logging.info(f"IUPAC\tIUPAC PROBLEM: {name} ({numberofIUPACs_1}, {numberofIUPACs_2})")
                ex.write(f"{name}\tIUPAC\tThe number of ambiguous bases ({numberofIUPACs_1}, {numberofIUPACs_2}) in sequence exceeds the number of allowed ambiguous bases (16, {allowed_number})\n")
        else:
            o.write(f">{name}\n")
            o.write(f"{seq}\n")
            o2.write(f">{name}\n")
            o2.write(f"{orig_seqs_dict[name]}\n")
