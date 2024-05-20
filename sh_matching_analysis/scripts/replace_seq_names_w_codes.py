import argparse
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to replace original seq. ids with unique internal ids.")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
infile = user_dir / f"source_{run_id}"
fasta_file = user_dir / f"source_{run_id}_fasta"
names_file = user_dir / f"source_{run_id}_names"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

seq_ct = 0
with open(infile, "r") as handle, open(fasta_file, "w") as f, open(names_file, "w") as n:
    for counter, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
        seq_ct += 1
        new_id = f"{run_id}_{counter}"

        f.write(f">i{new_id}i\n")
        f.write(f"{record.seq}\n")

        n.write(f"{record.id}\ti{new_id}i\t1\n")

logging.info(f"START\tNumber of sequences to start with: {seq_ct}")
