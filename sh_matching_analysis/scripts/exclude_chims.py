#!/usr/bin/python
import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to exclude chimeras identified by vsearch")
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
global_infile = user_dir / "usearch_global.full.75.blast6out.txt"
infile = user_dir / "seqs_out.fasta"
outfile = user_dir / "seqs_out_chim.fasta"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

nohit_counter = 0
nogo_counter = 0
chim_counter = 0
alveolates_counter = 0

len_limit_1 = 400
len_limit_2 = 350
if region == "its2":
    len_limit_1 = 100
    len_limit_2 = 100

global_chim_dict = {}

# open excluded seq file
with open(ex_file, "a") as ex, open(global_infile) as glob:
    # open blast6out file for parsing
    # TODO - csv.DictReader possibility
    dataReader_glob = csv.reader(glob, delimiter="\t")
    for row in dataReader_glob:
        # 1     Query label.
        # 2     Target (database sequence or cluster centroid) label.
        # 3     Percent identity.
        # 4     Alignment length.
        # 5     Number of mismatches.
        # 6     Number of gap opens.
        # 7     Start position in query. Query coordinates start with 1 at the first base in the sequence as it appears in the input file. For translated searches (nucleotide queries, protein targets), query start<end for +ve frame and start>end for -ve frame.
        # 8     End position in query.
        # 9     Start position in target. Target coordinates start with 1 at the first base in sequence as it appears in the database. For untranslated nucleotide searches, target start<end for plus strand, start>end for a reverse-complement alignment.
        # 10        End position in target.
        # 11        E-value calculated using Karlin-Altschul statistics.
        # 12        Bit score calculated using Karlin-Altschul statistics.
        if not row[1] == "*":
            len_of_query_cov = int(row[7]) - int(row[6])
            len_of_target_cov = int(row[9]) - int(row[8])
            len_alignment = int(row[3])

            perc_alignment_cov = len_alignment * 100 / len_of_query_cov
            perc_alignment_cov_target = len_alignment * 100 / len_of_target_cov

            if perc_alignment_cov <= 80 and perc_alignment_cov_target <= 80:
                global_chim_dict[row[0]] = 1
            elif perc_alignment_cov <= 85 and perc_alignment_cov_target <= 85:
                global_chim_dict[row[0]] = 1
            elif int(row[3]) < len_limit_1 and int(row[7]) >= len_limit_1:
                global_chim_dict[row[0]] = 1
        else:
            logging.info(f"CHIM\tNOHIT_RECORD\t{row[0]}")
            nohit_counter += 1

    logging.info(f"CHIM\tNo. of nohits: {nohit_counter}")

    with open(outfile, "w") as o, open(infile, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in global_chim_dict:
                if len(str(record.seq)) >= len_limit_2:
                    o.write(f">{record.id}\n")
                    o.write(f"{record.seq}\n")
                else:
                    nogo_counter += 1
                    logging.info(f"CHIM\tNOGO_RECORD\t{record.id}")
                    ex.write(f"{record.id}\tCHIM\tSequence too short (<{len_limit_2} nucl).\n")
            else:
                chim_counter += 1
                logging.info(f"CHIM\tCHIM_RECORD\t{record.id}")
                ex.write(f"{record.id}\tCHIM\tIdentified as chimeric based on the results of vsearch (usearch_global) chimera detection.\n")

    logging.info(f"CHIM\tNo go for {nogo_counter} (length) + {chim_counter} (chims) sequences.")
