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
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
itsx_dir = user_dir / "ITSx"
infile = itsx_dir / "itsx_sh_out.ITS2.full_and_partial.fasta"
pos_file = itsx_dir / "itsx_sh_out.positions.txt"
full_file = itsx_dir / "itsx_sh_out.full_and_partial.fasta"
outfile = user_dir / "seqs_out_1.fasta"

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
full_new_counter = 0  # TODO - doesn't seem to be used
its2_counter = 0

# read in ITS positions (to make sure that ITS1, 5.8S and ITS2 regions are all found, but may just be too short)
# open excluded seq file
with open(ex_file, "a") as ex, open(pos_file) as pos:
    # TODO - csv.DictReader would be nice
    dataReader_pos = csv.reader(pos, delimiter="\t")
    for row in dataReader_pos:
        if not row[3] == "ITS1: Not found" and not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
            chim_match = re.search("Chimeric!", row[7])
            if not chim_match:
                positions_dict[row[0]] = 1
                len_fields = row[1].split(" ")
                length_dict[row[0]] = int(len_fields[0])
            else:
                ex.write(f"{row[0]}\tPRINT_FAS\tChimeric according to ITSx.\n")

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
            if length_dict[record_id] >= 140:
                if record_id in full_seq_dict and len(full_seq_dict[record_id]) >= 140:
                    o.write(f">{record_id}\n")
                    o.write(f"{full_seq_dict[record_id]}\n")
                    full_counter += 1
                else:
                    if len(str(record.seq)) >= 140:
                        o.write(f">{record_id}\n")
                        o.write(f"{record.seq}\n")
                        its2_counter += 1
            else:
                logging.info(f"PRINT_FAS\tSequence too short - {record_id}")
                ex.write(f"{record_id}\tPRINT_FAS\tSequence too short.\n")
        else:
            no_coverage_count += 1

logging.info(f"PRINT_FAS\tNo coverage for {no_coverage_count} sequences.")
logging.info(
    f"PRINT_FAS\tNo. of full seqs: {full_counter}; "
    f"No. of full new seqs: {full_new_counter}; "
    f"No of ITS2 seqs: {its2_counter}"
)
