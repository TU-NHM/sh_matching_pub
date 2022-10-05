import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to create mapping file for duplicates")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")

# infiles
tmp_file = user_dir / "clusters_pre" / "tmp.txt"
singleton_file = user_dir / "clusters_pre" / "singletons.txt"
cov100_uniq_file = user_dir / f"source_{run_id}_fastanames"
cov96_uniq_file = user_dir / "clusters_100.uc"

# outfiles
duplic_seqs_file = user_dir / "duplic_seqs.txt"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

seq_counter = 0  # sequence count
seq_counter_all = 0  # sequence count with duplicates included

# include duplicate sequences (tmp_files/sequences.names)
print("Collecting duplicate sequences (cov100) ...")
cov100_uniq_dict = {}
with open(cov100_uniq_file, "r") as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        # include only those sequences where duplicates are present
        if row[0] != row[1]:
            cov100_list = row[1].split(",")
            cov100_list.remove(row[0])
            cov100_uniq_dict[row[0]] = (",").join(cov100_list)

# include vsearch 100%sim/96%cov cluster members (clusters_100.uc)
print("Collecting duplicate sequences (cov96) ...")
cov96_uniq_dict = {}
with open(cov96_uniq_file, "r") as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        if row[0] == "H":
            if row[9] in cov96_uniq_dict:
                if row[8] in cov100_uniq_dict:
                    cov96_uniq_dict[row[9]] = cov96_uniq_dict[row[9]] + "," + row[8] + ","+ cov100_uniq_dict[row[8]]
                else:
                    cov96_uniq_dict[row[9]] = cov96_uniq_dict[row[9]] + "," + row[8]
            else:
                if row[8] in cov100_uniq_dict and row[9] in cov100_uniq_dict:
                    cov96_uniq_dict[row[9]] = cov100_uniq_dict[row[9]] + "," + row[8] + "," + cov100_uniq_dict[row[8]]
                elif row[9] in cov100_uniq_dict:
                    cov96_uniq_dict[row[9]] = cov100_uniq_dict[row[9]] + "," + row[8]
                elif row[8] in cov100_uniq_dict:
                    cov96_uniq_dict[row[9]] = row[8] + "," + cov100_uniq_dict[row[8]]
                else:
                    cov96_uniq_dict[row[9]] = row[8]
        elif row[0] == "S":
            if not row[8] in cov96_uniq_dict:
                if row[8] in cov100_uniq_dict:
                    cov96_uniq_dict[row[8]] = cov100_uniq_dict[row[8]]

with open(duplic_seqs_file, "w") as dupl:
    # process alignments folder
    with open(tmp_file, "r") as t:
        dataReader = csv.reader(t, delimiter="\t")
        for row in dataReader:
            name = row[0]
            ex_cl_file = user_dir / "clusters_pre" / "clusters" / name
            print("Processing " + name)
            with open(ex_cl_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq_counter += 1
                    seq_counter_all += 1
                    seq_list = []
                    if record.id in cov96_uniq_dict:
                        seq_list = cov96_uniq_dict[record.id].split(",")
                    for seq in seq_list:
                        if not seq == record.id:
                            seq_counter_all += 1
                            parent_seq = str(record.id)
                            dupl.write(str(seq) + "\t" + str(parent_seq) + "\t" + name + "\n")

            print("Finished with " + name)
            print("\n\nNumber of sequences at the end: " + str(seq_counter) + "\n\n")


    # process singletons folder
    with open(singleton_file, "r") as s:
        dataReader = csv.reader(s, delimiter="\t")
        for row in dataReader:
            name = row[0]
            print("Processing " + name)
            ex_singleton_file = user_dir / "clusters_pre" / "singletons" / name
            with open(ex_singleton_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq_counter += 1
                    seq_counter_all += 1
                    seq_list = []
                    if record.id in cov96_uniq_dict:
                        seq_list = cov96_uniq_dict[record.id].split(",")
                    for seq in seq_list:
                        if not seq == record.id:
                            seq_counter_all += 1
                            parent_seq = str(record.id)
                            dupl.write(str(seq) + "\t" + str(parent_seq) + "\t" + name + "\n")

            print("Finished with " + name)
            print("\n\nNumber of sequences at the end: " + str(seq_counter) + "\n\n")

logging.info(f"USEARCH_PARSER\tNumber of sequences at the end: {seq_counter} and {seq_counter_all}")
