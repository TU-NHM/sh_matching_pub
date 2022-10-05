#!/usr/bin/python
import argparse
import logging
import os
import re
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to select representative sequences from vsearch 100 percent clustering")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile_centroids = user_dir / "centroids_100.fasta"
infile_iupac = user_dir / "iupac_out_full.fasta"
outfile = user_dir / "iupac_out_vsearch.fasta"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)
ex_file = user_dir / f"excluded_{run_id}.txt"

seq_counter = 0
disc_seq_counter = 0

# read in all sequences that came out from ITSx & IUPAC QC steps
iupac_full_seqs_dict = {}
with open(infile_iupac, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        iupac_full_seqs_dict[record.id] = record.seq

# go through vsearch centroids file to create a new dataset of vsearch representatives
with open(outfile, "w") as o, open(infile_centroids, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in iupac_full_seqs_dict:
            seq_counter += 1
            o.write(">" + str(record.id) + "\n" + str(iupac_full_seqs_dict[record.id]) + "\n")
        else:
            ex.write(f"{record.id}\tVSEARCH\tDiscarded during the VSEARCH 100% similarity/96% coverage clustering step")
            disc_seq_counter += 1

logging.info(f"Number of sequences passed vsearch RepS selection: {seq_counter}")
logging.info(f"Number of sequences discarded: {disc_seq_counter}")
