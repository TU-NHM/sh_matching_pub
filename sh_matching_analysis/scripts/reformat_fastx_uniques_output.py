#!/usr/bin/python
import argparse
import csv
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format vsearch --fastx_uniques output to restore original sequence.")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile = user_dir / f"source_{run_id}_fastauc"
outfile = user_dir / f"source_{run_id}_fastanames"

matches_dict = {}
with open(infile, "r") as i:
    dataReader = csv.reader(i, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            matches_dict[row[8]] = row[8]
        elif row[0] == "H":
            matches_dict[row[9]] = matches_dict[row[9]] + "," + row[8]

with open(outfile, "w") as o:
    for key in matches_dict:
        o.write(f"{key}\t{matches_dict[key]}\n")
