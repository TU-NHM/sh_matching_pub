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
no_detect_file_o = itsx_dir_o / "itsx_sh_out_o_no_detections.txt"
no_detect_file_f = itsx_dir_f / "itsx_sh_out_no_detections.txt"
full_file = itsx_dir_o / "itsx_sh_out_o.full_and_partial.fasta"
outfile = user_dir / "seqs_out_2.fasta"

cov100_uniq_file = user_dir / f"source_{run_id}_fastanames"

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
its2_counter = 0
len_limit = 140
ex_o_ct = 0

# include info about duplicate sequences
cov100_uniq_dict = {}
with open(cov100_uniq_file, "r") as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        # include only those sequences where duplicates are present
        cov100_uniq_dict[row[0]] = row[1]

# read in ITS positions (to make sure that ITS1, 5.8S and ITS2 regions are all found, but may just be too short)
if region == "itsfull":
    with open(pos_file_f) as pos:
        # TODO - csv.DictReader possibility
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
elif region == "its2":
    len_limit = 100
    with open(pos_file_f) as pos:
        dataReader_pos = csv.reader(pos, delimiter="\t")
        for row in dataReader_pos:
            if not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not chim_match and not too_long_match:
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
    with open(pos_file_o) as pos_o:
        # TODO - csv.DictReader possibility
        dataReader_pos_o = csv.reader(pos_o, delimiter="\t")
        for row in dataReader_pos_o:
            if not row[3] == "ITS1: Not found" and not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                part_58S_match = re.search("Broken or partial sequence, only partial 5.8S!", row[7])
                no_58S_match = re.search("Broken or partial sequence, no 5.8S!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not chim_match and not part_58S_match and not no_58S_match and not too_long_match:
                    new_positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    new_length_dict[row[0]] = int(len_fields[0])
                else:
                    logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS_O\tChimeric or broken sequence according to ITSx.")
            else:
                ex_o_ct += 1
                logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS_O\tITS1 or ITS2 sequence not detected.")
elif region == "its2":
    with open(pos_file_o) as pos_o:
        # TODO - csv.DictReader possibility
        dataReader_pos_o = csv.reader(pos_o, delimiter="\t")
        for row in dataReader_pos_o:
            if not row[5] == "ITS2: Not found" and not row[5] == "ITS2: No start" and not row[5] == "ITS2: No end":
                chim_match = re.search("Chimeric!", row[7])
                too_long_match = re.search("ITS region too long!", row[7])
                if not chim_match and not too_long_match:
                    new_positions_dict[row[0]] = 1
                    len_fields = row[1].split(" ")
                    new_length_dict[row[0]] = int(len_fields[0])
                else:
                    logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS_O\tChimeric or broken sequence according to ITSx.")
            else:
                ex_o_ct += 1
                logging.info(f"{cov100_uniq_dict[row[0]]}\tPRINT_FAS_O\tITS1 or ITS2 sequence not detected.")

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
                    logging.info(f"PRINT_FAS_O\tSequence too short - {cov100_uniq_dict[record_id]}")
                    ex.write(f"{cov100_uniq_dict[record_id]}\tPRINT_FAS_O\tSequence too short.\n")
            else:
                no_coverage_count += 1
                ex.write(f"{cov100_uniq_dict[record_id]}\tPRINT_FAS_O\tITS1 or ITS2 sequence not detected.\n")

logging.info(f"PRINT_FAS_O\tNo coverage for {no_coverage_count} sequences.")
logging.info(
    f"PRINT_FAS_O\tNo. of full seqs: {full_counter}; "
    f"No of ITS2 seqs: {its2_counter}"
)

# read in no_detections.txt files for both runs and print out info about excluded sequences
no_detect_f_set = set()
if no_detect_file_f.is_file():
    with open (no_detect_file_f, "r") as det:
        dataReader = csv.reader(det, delimiter="\t")
        for row in dataReader:
            no_detect_f_set.add(row[0])

if no_detect_file_o.is_file():
    with open(ex_file, "a") as ex, open (no_detect_file_o, "r") as det:
        dataReader = csv.reader(det, delimiter="\t")
        for row in dataReader:
            if row[0] in no_detect_f_set:
                ex.write(f"{cov100_uniq_dict[record_id]}\tPRINT_FAS_O\tSequence not detected as ITS by ITSx.\n")

logging.info(
    f"No of seqs (other) excluded (no ITS1 or ITS2 region detected): {ex_o_ct}"
)
