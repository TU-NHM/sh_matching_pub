import argparse
import csv
import logging
import os
import sys
from pathlib import Path

from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description="Script to parse USEARCH clustering output (90) for next step SH clustering")
parser.add_argument("run_id", help="Need run id in numeric format!")
parser.add_argument("name", help="Need cluster name!")
args = parser.parse_args()

# read in args
run_id = args.run_id
name = args.name
name_folder = name + "_folder"
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
file = user_dir / "compounds" / "clusters_2_90.uc"
tmp_file1 = user_dir / "clusters_out_2_90.txt"
tmp_file_nohits = user_dir / "iupac_out_vsearch.fasta"
tmp_file_compound = user_dir / "compounds" / name
tmp_cl_file = user_dir / "compounds" / name_folder / "tmp.txt"
tmp_singl_file = user_dir / "compounds" / name_folder / "singletons.txt"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

cluster_dict = {}
cluster_counter = 0
seq_counter = 0
original_seq_dict = {}

with open(file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            # seed sequence, create new cluster
            cluster_dict[row[1]] = f"{cluster_counter}\t{row[1]}\t"
            cluster_dict[row[1]] = cluster_dict[row[1]] + row[8] + " "
            seq_counter += 1
            cluster_counter += 1
        elif row[0] == "H":
            # hit with target sequence
            # print previous rounds
            cluster_dict[row[1]] = cluster_dict[row[1]] + row[8] + " "
            seq_counter += 1
        elif row[0] == "C":
            # cluster centroid, ignore at the moment (same as H)
            continue
        else:
            logging.info(f"CLP_F_90\t{row[1]}")

logging.info(f"CLP_F_90\tTotal no of sequences: {seq_counter}")

with open(tmp_file1, "w") as o:
    for key, value in cluster_dict.items():
        cluster_dict[key] = value.strip()
        o.write(cluster_dict[key] + "\n")

# create hash for sequences
with open(tmp_file_compound, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_seq_dict[record.id] = str(record.seq)
# add users' sequences
with open(tmp_file_nohits, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_seq_dict[record.id] = str(record.seq)

# start dividing clusters into real clusters and singletons in fasta files
with open (tmp_cl_file, "w") as cl, open(tmp_singl_file, "w") as singl, open(tmp_file1) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        cluster_seqs = row[2].split(" ")
        if len(cluster_seqs) == 1:
            # singleton cluster
            singl_file = user_dir / "compounds" / name_folder / "singletons" / f"Singleton{row[0]}"
            with open(singl_file, "w") as s:
                s.write(f">{cluster_seqs[0]}\n")
                s.write(f"{original_seq_dict[cluster_seqs[0]]}\n")
            singl.write(f"Singleton{row[0]}\n")
        else:
            cl_file = user_dir / "compounds" / name_folder / "clusters" / f"Cluster{row[0]}"
            with open(cl_file, "w") as c:
                for item in cluster_seqs:
                    c.write(f">{item}\n{original_seq_dict[item]}\n")
            cl.write(f"Cluster{row[0]}\n")
