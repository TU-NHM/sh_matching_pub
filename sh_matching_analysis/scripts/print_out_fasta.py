#!/usr/bin/python
import argparse
import csv
import logging
import os
import re
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format ITSx output (fungi)")
parser.add_argument("run_id", help="Need run id in numeric format!")
parser.add_argument("region", help="Need region (either its2 or itsfull)!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

region = args.region
if not region == "its2" and not region == "itsfull":
    raise ValueError("Region is not one of: its2, itsfull", region)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
itsx_dir = user_dir / "ITSx"
infile = itsx_dir / "itsx_sh_out.ITS2.full_and_partial.fasta"
pos_file = itsx_dir / "itsx_sh_out.positions.txt"
full_file = itsx_dir / "itsx_sh_out.full_and_partial.fasta"
outfile = user_dir / "seqs_out_1.fasta"

cov100_uniq_file = user_dir / f"source_{run_id}_fastanames"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

positions_dict = {}
length_dict = {}
no_coverage_count = 0
full_seq_dict = {}
full_counter = 0
its2_counter = 0
len_limit = 140
ex_f_ct = 0

# include info about duplicate sequences
cov100_uniq_dict = {}
with open(cov100_uniq_file, "r") as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        # include only those sequences where duplicates are present
        cov100_uniq_dict[row[0]] = row[1]

# read in ITS positions (to make sure that ITS1, 5.8S and ITS2 regions are all found, but may just be too short)
# open excluded seq file
if region == "itsfull":
    with open(ex_file, "a") as ex, open(pos_file) as pos:
        # TODO - csv.DictReader would be nice
        dataReader_pos = csv.reader(pos, delimiter="\t")
        for row in dataReader_pos:
            if not row[3] == "ITS1: Not found" and not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                part_58S_match = re.search("Broken or partial sequence, only partial 5.8S!", row[7])
                no_58S_match = re.search("Broken or partial sequence, no 5.8S!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not chim_match and not part_58S_match and not no_58S_match and not too_long_match:
                    positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    length_dict[row[0]] = int(len_fields[0])
                else:
                    ex.write(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS\tChimeric or broken sequence according to ITSx.\n")
            else:
                ex_f_ct += 1
                logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS\tITS1 or ITS2 sequence not detected.")
elif region == "its2":
    len_limit = 100
    with open(ex_file, "a") as ex, open(pos_file) as pos:
        # TODO - csv.DictReader would be nice
        dataReader_pos = csv.reader(pos, delimiter="\t")
        for row in dataReader_pos:
            if not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not chim_match and not too_long_match:
                    positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    length_dict[row[0]] = int(len_fields[0])
                else:
                    ex.write(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS\tChimeric or broken sequence according to ITSx.\n")
            else:
                ex_f_ct += 1
                logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS\tITS1 or ITS2 sequence not detected.\n")

with open(full_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        name_fields = name.split("|")
        full_seq_dict[name_fields[0]] = str(record.seq)

# open excluded seq file
with open(ex_file, "a") as ex, open(outfile, "w") as o, open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        name_fields = name.split("|")
        record_id = name_fields[0]
        if record_id in positions_dict:
            if length_dict[record_id] >= len_limit:
                if record_id in full_seq_dict and len(full_seq_dict[record_id]) >= len_limit:
                    o.write(f">{record_id}\n")
                    o.write(f"{full_seq_dict[record_id]}\n")
                    full_counter += 1
                else:
                    if len(str(record.seq)) >= len_limit:
                        o.write(f">{record_id}\n")
                        o.write(f"{record.seq}\n")
                        its2_counter += 1
            else:
                logging.info(f"PRINT_FAS\tSequence too short - {cov100_uniq_dict[record_id]}")
                ex.write(f"{cov100_uniq_dict[record_id]}\tPRINT_FAS\tSequence too short.\n")
        else:
            no_coverage_count += 1

logging.info(f"PRINT_FAS\tNo coverage for {no_coverage_count} sequences.")
logging.info(
    f"PRINT_FAS\tNo. of full seqs: {full_counter}; "
    f"No of ITS2 seqs: {its2_counter}"
)
logging.info(
    f"No of seqs (fungi) excluded (no ITS region detected): {ex_f_ct}"
)
