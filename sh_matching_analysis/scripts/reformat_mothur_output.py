#!/usr/bin/python
import argparse
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format mothur output to restore original sequence.")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile = user_dir / f"source_{run_id}_fasta"
mothur = user_dir / f"source_{run_id}_fastaunique_mothur"
mothur_orig = user_dir / f"source_{run_id}_fastaunique"

original_dict = {}
with open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_dict[record.id] = str(record.seq)

with open(mothur, "r") as handle, open(mothur_orig, "w") as f:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in original_dict:
            f.write(f">{record.id}\n")
            f.write(f"{record.seq}\n")
