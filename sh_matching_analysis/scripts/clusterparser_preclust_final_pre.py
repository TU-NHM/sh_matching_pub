import argparse
import csv
import logging
import os
import sys
from pathlib import Path

from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description="Script to parse USEARCH preclustering outputs (all)")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
file = user_dir / "clusters_80_pre.uc"
tmp_file1 = user_dir / "clusters_out_97_pre.txt"
tmp_file2 = user_dir / "clusters_out_95_pre.txt"
tmp_file3 = user_dir / "clusters_out_90_pre.txt"
tmp_file4 = user_dir / "clusters_out_80_pre.txt"
tmp_file_nohits = user_dir / "iupac_out_vsearch.fasta"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

round1_dict = {}
round2_dict = {}
round3_dict = {}
cluster_dict = {}
cluster_counter = 0
seq_counter = 0
original_seq_dict = {}

# get round1 clusters
with open(tmp_file1) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        round1_dict[row[2]] = row[3]

# get round2 clusters
with open(tmp_file2) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        round2_dict[row[2]] = row[3]

# get round3 clusters
with open(tmp_file3) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        round3_dict[row[2]] = row[3]

with open(file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        if row[0] == "S":
            # seed sequence, create new cluster
            cluster_dict[row[1]] = f"{cluster_counter}\t{row[1]}\t"
            # print previous rounds
            if row[8] in round3_dict:
                tmp_fields3 = round3_dict[row[8]].split(" ")
                for item3 in tmp_fields3:
                    if item3 in round2_dict:
                        tmp_fields2 = round2_dict[item3].split(" ")
                        for item2 in tmp_fields2:
                            if item2 in round1_dict:
                                tmp_fields1 = round1_dict[item2].split(" ")
                                for item1 in tmp_fields1:
                                    cluster_dict[row[1]] = cluster_dict[row[1]] + item1 + " "
                                    seq_counter += 1
                            else:
                                cluster_dict[row[1]] = cluster_dict[row[1]] + item2 + " "
                                seq_counter += 1
                    else:
                        cluster_dict[row[1]] = cluster_dict[row[1]] + item3 + " "
                        seq_counter += 1
            else:
                cluster_dict[row[1]] = cluster_dict[row[1]] + row[8] + " "
                seq_counter += 1
            cluster_counter += 1
        elif row[0] == "H":
            # hit with target sequence
            # print previous rounds
            if row[8] in round3_dict:
                tmp_fields3 = round3_dict[row[8]].split(" ")
                for item3 in tmp_fields3:
                    if item3 in round2_dict:
                        tmp_fields2 = round2_dict[item3].split(" ")
                        for item2 in tmp_fields2:
                            if item2 in round1_dict:
                                tmp_fields1 = round1_dict[item2].split(" ")
                                for item1 in tmp_fields1:
                                    cluster_dict[row[1]] = cluster_dict[row[1]] + item1 + " "
                                    seq_counter += 1
                            else:
                                cluster_dict[row[1]] = cluster_dict[row[1]] + item2 + " "
                                seq_counter += 1
                    else:
                        cluster_dict[row[1]] = cluster_dict[row[1]] + item3 + " "
                        seq_counter += 1
            else:
                cluster_dict[row[1]] = cluster_dict[row[1]] + row[8] + " "
                seq_counter += 1
        elif row[0] == "C":
            # cluster centroid, ignore at the moment (same as H)
            continue
        else:
            logging.info(f"CLPP_F\t{row[1]}")

logging.info(f"CLPP_F\tTotal no of sequences: {seq_counter}")

with open(tmp_file4, "w") as o:
    for key, value in cluster_dict.items():
        cluster_dict[key] = value.strip()
        o.write(cluster_dict[key] + "\n")

# create dict for sequences
with open(tmp_file_nohits, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        original_seq_dict[record.id] = str(record.seq)

# start dividing clusters into real clusters and singletons in fasta files
with open(tmp_file4, "r") as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        cluster_seqs = row[2].split(" ")
        if len(cluster_seqs) == 1:
            # singleton cluster
            singl_file = user_dir / "clusters_pre" / "singletons" / f"Singleton{row[0]}"
            with open(singl_file, "w") as s:
                s.write(f">{cluster_seqs[0]}\n")
                s.write(f"{original_seq_dict[cluster_seqs[0]]}\n")
        else:
            cl_file = user_dir / "clusters_pre" / "clusters" / f"Cluster{row[0]}"
            with open(cl_file, "w") as c:
                for item in cluster_seqs:
                    c.write(f">{item}\n{original_seq_dict[item]}\n")
