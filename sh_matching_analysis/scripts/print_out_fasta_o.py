#!/usr/bin/python
import argparse
import csv
import logging
import os
import re
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to format ITSx output (others)")
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
itsx_dir_o = user_dir / "ITSx_o"
itsx_dir_f = user_dir / "ITSx"
infile_o = itsx_dir_o / "itsx_sh_out_o.ITS2.full_and_partial.fasta"
infile_f = itsx_dir_f / "itsx_sh_out.ITS2.full_and_partial.fasta"
pos_file_o = itsx_dir_o / "itsx_sh_out_o.positions.txt"
pos_file_f = itsx_dir_f / "itsx_sh_out.positions.txt"
full_file = itsx_dir_o / "itsx_sh_out_o.full_and_partial.fasta"
outfile = user_dir / "seqs_out_2.fasta"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

fungi_dict = {}
positions_dict = {}
length_dict = {}
new_positions_dict = {}
new_length_dict = {}
no_coverage_count = 0
full_seq_dict = {}
full_counter = 0
full_new_counter = 0  # TODO - doesn't seem to be used
its2_counter = 0
len_limit = 140

# read in ITS positions (to make sure that ITS1, 5.8S and ITS2 regions are all found, but may just be too short)
if region == "itsfull":
    with open(pos_file_f) as pos:
        # TODO - csv.DictReader possibility
        dataReader_pos = csv.reader(pos, delimiter="\t")
        for row in dataReader_pos:
            if not row[3] == "ITS1: Not found" and not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                if not chim_match:
                    positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    length_dict[row[0]] = int(len_fields[0])
elif region == "its2":
    len_limit = 100
    with open(pos_file_f) as pos:
        # TODO - csv.DictReader possibility
        dataReader_pos = csv.reader(pos, delimiter="\t")
        for row in dataReader_pos:
            if not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                if not chim_match:
                    positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    length_dict[row[0]] = int(len_fields[0])

with open(infile_f, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        name_fields = name.split("|")
        fungi_dict[name_fields[0]] = 1

# open excluded seq file
if region == "itsfull":
    with open(ex_file, "a") as ex, open(pos_file_o) as pos_o:
        # TODO - csv.DictReader possibility
        dataReader_pos_o = csv.reader(pos_o, delimiter="\t")
        for row in dataReader_pos_o:
            if not row[3] == "ITS1: Not found" and not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                if not chim_match:
                    new_positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    new_length_dict[row[0]] = int(len_fields[0])
                else:
                    ex.write(f"{row[0]}\tPRINT_FAS_O\tChimeric according to ITSx.\n")
elif region == "its2":
    with open(ex_file, "a") as ex, open(pos_file_o) as pos_o:
        # TODO - csv.DictReader possibility
        dataReader_pos_o = csv.reader(pos_o, delimiter="\t")
        for row in dataReader_pos_o:
            if not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                if not chim_match:
                    new_positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    new_length_dict[row[0]] = int(len_fields[0])
                else:
                    ex.write(f"{row[0]}\tPRINT_FAS_O\tChimeric according to ITSx.\n")

# get ITS sequence from full file into hash
with open(full_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        name_fields = name.split("|")
        full_seq_dict[name_fields[0]] = str(record.seq)

# open excluded seq file
with open(ex_file, "a") as ex, open(outfile, "w") as o, open(infile_o, "r") as handle:
    # print out fasta file with correct headers
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        name_fields = name.split("|")
        record_id = name_fields[0]
        if record_id not in fungi_dict:
            if record_id in new_positions_dict:
                if new_length_dict[record_id] >= len_limit:
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
                    logging.info(f"PRINT_FAS_O\tSequence too short - {record_id}")
                    ex.write(f"{record_id}\tPRINT_FAS_O\tSequence too short.\n")
            else:
                no_coverage_count += 1

logging.info(f"PRINT_FAS_O\tNo coverage for {no_coverage_count} sequences.")
logging.info(
    f"PRINT_FAS_O\tNo. of full seqs: {full_counter}; "
    f"No. of full new seqs: {full_new_counter}; "
    f"No of ITS2 seqs: {its2_counter}"
)
