import argparse
import csv
import glob
import os
import re
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to parse blastclust output")
parser.add_argument("run_id", help="Need run id in numeric format!")
parser.add_argument("threshold", help="Need threshold numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
threshold = args.threshold
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)
if not threshold.isdigit():
    raise ValueError("Threshold is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
matches_dir = user_dir / "matches"
matches_file = matches_dir / f"matches_1_{threshold}.txt"

prev_dict = {}

# read in results from prev version (if singleton, no need to check it again)
prev_file = None
if threshold == "975":
    prev_file = matches_dir / "matches_1_97.txt"
elif threshold == "98":
    prev_file = matches_dir / "matches_1_975.txt"
elif threshold == "985":
    prev_file = matches_dir / "matches_1_98.txt"
elif threshold == "99":
    prev_file = matches_dir / "matches_1_985.txt"
elif threshold == "995":
    prev_file = matches_dir / "matches_1_99.txt"
elif threshold == "100":
    prev_file = matches_dir / "matches_1_995.txt"

if prev_file:
    with open(prev_file) as p:
        dataReader = csv.reader(p, delimiter="\t")
        for row in dataReader:
            if row[2] == "singleton":
                prev_dict[row[0]] = 1

# open matches file
with open(matches_file, "w") as o:
    # go through BC output
    glob_match = f"{user_dir}/blastclust_1/Cluster*"
    file_list = glob.glob(glob_match)
    for file in file_list:
        if re.search(r'Cluster\d{1,}$', file) is not None:
            tmp_bc_file = f"{file}_out_{threshold}"
            ucl_code = file.split("/")[-1]
            with open(tmp_bc_file) as f:
                dataReader = csv.reader(f, delimiter=" ")
                for row in dataReader:
                    if row[0] not in prev_dict:
                        if len(row) == 1:
                            o.write(f"{row[0]}\t\tsingleton\t{ucl_code}\t\n")
                        else:
                            joined_row = " ".join(row)
                            for seq in row:
                                if not seq == "":
                                    o.write(f"{seq}\t\tnew cluster\t{ucl_code}\t{joined_row}\n")
                    else:
                        o.write(f"{row[0]}\t\tsingleton\t{ucl_code}\t\n")

    # process singletons folder
    glob_match_s = f"{user_dir}/clusters/singletons/Singleton*"
    file_list_s = glob.glob(glob_match_s)
    for file_s in file_list_s:
        ucl_code = file_s.split("/")[-1]
        with open(file_s, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id not in prev_dict:
                    o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")
                else:
                    o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")
