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
cov100_uniq_file = user_dir / f"source_{run_id}_fastanames"

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
with open(cov100_uniq_file, "r") as f, open(duplic_seqs_file, "w") as dupl:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        seq_counter += 1
        # include only those sequences where duplicates are present
        if row[0] != row[1]:
            cov100_list = row[1].split(",")
            cov100_list.remove(row[0])
            for seq in cov100_list:
                seq_counter_all += 1
                dupl.write(str(seq) + "\t" + str(row[0]) + "\t\n")

logging.info(f"USEARCH_PARSER_UNIQ\tNumber of sequences at the end: {seq_counter} and {seq_counter_all}")
