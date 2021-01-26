import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to remove sequences with more than X ambiguous bases)")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile = user_dir / "iupac_out.fasta"
outfile1 = user_dir / "hits.fasta"
outfile2 = user_dir / "nohits.fasta"
outfile3 = user_dir / "hits.txt"
map_file = user_dir / "closedref.75.map.uc"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

full_dict = {}
with open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        full_dict[name] = str(record.seq)

with open(outfile1, "w") as o1, open(outfile2, "w") as o2, open(outfile3, "w") as o3, open(map_file) as map_f:
    # TODO - possibility for csv.DictReader
    dataReader_map = csv.reader(map_f, delimiter="\t")
    for row in dataReader_map:
        if row[0] == "N":
            o2.write(f">{row[8]}\n")
            o2.write(f"{full_dict[row[8]]}\n")
        elif row[0] == "H":
            o3.write(f"{row[8]}\t{row[2]}\t{row[3]}\t{row[9]}\n")
            o1.write(f">{row[8]}\n")
            o1.write(f"{full_dict[row[8]]}\n")
        else:
            logging.info(f"VSEARCH\tERR: {row[0]}")
