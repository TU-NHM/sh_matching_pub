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

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

matches_dict = {}
with open(infile, "r") as i:
    dataReader = csv.reader(i, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            matches_dict[row[8]] = row[8]
        elif row[0] == "H":
            matches_dict[row[9]] = matches_dict[row[9]] + "," + row[8]

uniq_ct = 0
with open(outfile, "w") as o:
    for key in matches_dict:
        uniq_ct += 1
        o.write(f"{key}\t{matches_dict[key]}\n")

logging.info(f"FASTX_UNIQUES\tNumber of unique sequences: {uniq_ct}")
