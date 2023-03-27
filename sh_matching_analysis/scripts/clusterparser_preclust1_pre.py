import argparse
import csv
import logging
import os
import sys
from pathlib import Path

from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description="Script to parse USEARCH preclustering output (97) for next clustering step (95)")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
file = user_dir / "clusters_97_pre.uc"
tmp_file1 = user_dir / "clusters_out_97_pre.txt"
tmp_file2 = user_dir / "iupac_out_vsearch.fasta"
tmp_file3 = user_dir / "in_95_pre.fasta"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

cluster_dict = {}
cluster_count_dict = {}
cluster_counter = 0
seq_counter = 0
original_seq_dict = {}

with open(file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            # seed sequence, create new cluster
            cluster_dict[row[1]] = f"{cluster_counter}\t{row[1]}\t{row[8]}\t{row[8]}"
            cluster_count_dict[row[1]] = 1
            cluster_counter += 1
            seq_counter += 1
        elif row[0] == "H":
            # hit with target sequence
            cluster_dict[row[1]] = cluster_dict[row[1]] + " " + row[8]
            cluster_count_dict[row[1]] += 1
            seq_counter += 1
        elif row[0] == "C":
            # cluster centroid, ignore at the moment (same as H)
            continue
        else:
            logging.info(f"CLP_1\t{row[1]}\n")

logging.info(f"CLPP_1\tTotal no of sequences: {seq_counter}")

# create dict for sequences
with open(tmp_file2, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_seq_dict[record.id] = str(record.seq)

# create input file for round2 clustering (95%)
with open(tmp_file1, "w") as o1, open(tmp_file3, "w") as o3:
    for key, value in cluster_dict.items():
        o1.write(f"{value}\n")
        row = value.split("\t")
        cluster_seqs = row[3].split(" ")
        o3.write(f">{cluster_seqs[0]}\n")
        o3.write(f"{original_seq_dict[cluster_seqs[0]]}\n")
