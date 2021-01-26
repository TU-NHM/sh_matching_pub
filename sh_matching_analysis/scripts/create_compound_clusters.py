import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to create compound clusters for hits")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
data_dir = Path("/sh_matching/data")
infile = user_dir / "iupac_out.fasta"
tmp_infile = user_dir / "closedref.75.map.uc"
compounds_file = user_dir / "compounds.txt"
exact_matches_file = user_dir / "exact_matches.txt"
sanger_refs_file = data_dir / "sanger_refs_sh.fasta"
compound2seq_file = data_dir / "compound2seq_mapping.txt"
sh2compound_file = data_dir / "sh2compound_mapping.txt"


# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)


exact_match_counter = 0
uniques = {}
sh_ucl_dict = {}
seq_ucl_dict = {}
seq_ucl_length_hash = {}
seq_ucl_map_hash = {}
original = {}
ucl_seq_dict = {}
seq_ucl_map_dict = {}
seq_ucl_length_dict = {}
must_refs_dict = {}

with open(sanger_refs_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        fields = name.split("_")
        new_name = fields[0]
        uniques[new_name] = str(record.seq)
        must_refs_dict[new_name] = 1

with open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        uniques[record.id] = str(record.seq)

# get compound and SH mappings
with open(sh2compound_file) as sh2c:
    dataReader_sh2c = csv.DictReader(sh2c, delimiter="\t", fieldnames=["field_1", "field_2", "field_3"])
    for row in dataReader_sh2c:
        sh_ucl_dict[row["field_1"]] = row["field_2"]

# get compound and sequence mappings - take "nearly identical" fact in account
with open(compound2seq_file) as c2seq:
    # TODO - DictReader possibility
    dataReader_c2seq = csv.reader(c2seq, delimiter="\t")
    for row in dataReader_c2seq:
        new_name = f"i{row[1]}i"

        if row[0] in seq_ucl_length_dict:
            if row[4] in seq_ucl_length_dict[row[0]]:
                if new_name in must_refs_dict:
                    if seq_ucl_map_dict[row[0]][row[4]] in must_refs_dict:
                        # add prev value to seq_ucl_dict
                        if row[0] in seq_ucl_dict:
                            seq_ucl_dict[row[0]] = f"{seq_ucl_dict[row[0]]},{seq_ucl_map_dict[row[0]][row[4]]}"
                        else:
                            seq_ucl_dict[row[0]] = seq_ucl_map_dict[row[0]][row[4]]
                    seq_ucl_map_dict[row[0]][row[4]] = new_name
                    seq_ucl_length_dict[row[0]][row[4]] = row[3]
                else:
                    if int(seq_ucl_length_dict[row[0]][row[4]]) < int(row[3]):
                        if seq_ucl_map_dict[row[0]][row[4]] in must_refs_dict:
                            # add prev value to seq_ucl_dict
                            if row[0] in seq_ucl_dict:
                                seq_ucl_dict[row[0]] = f"{seq_ucl_dict[row[0]]},{seq_ucl_map_dict[row[0]][row[4]]}"
                            else:
                                seq_ucl_dict[row[0]] = seq_ucl_map_dict[row[0]][row[4]]
                        seq_ucl_map_dict[row[0]][row[4]] = new_name
                        seq_ucl_length_dict[row[0]][row[4]] = row[3]
            else:
                seq_ucl_map_dict[row[0]][row[4]] = new_name
                seq_ucl_length_dict[row[0]][row[4]] = row[3]
        else:
            seq_ucl_map_dict[row[0]] = {}
            seq_ucl_map_dict[row[0]][row[4]] = new_name
            seq_ucl_length_dict[row[0]] = {}
            seq_ucl_length_dict[row[0]][row[4]] = row[3]

        original[new_name] = row[2]
        # ucl = row[0]
        # if ucl in seq_ucl_dict:
        #     seq_ucl_dict[ucl] = f"{seq_ucl_dict[ucl]},{new_name}"
        # else:
        #     seq_ucl_dict[ucl] = new_name

for ucl in seq_ucl_map_dict:
    for sh in seq_ucl_map_dict[ucl]:
        if ucl in seq_ucl_dict:
            seq_ucl_dict[ucl] = f"{seq_ucl_dict[ucl]},{seq_ucl_map_dict[ucl][sh]}"
        else:
            seq_ucl_dict[ucl] = seq_ucl_map_dict[ucl][sh]

with open(exact_matches_file, "w") as o, open(tmp_infile) as tmp_inf:
    dataReader_tmp_inf = csv.reader(tmp_inf, delimiter="\t")
    for row in dataReader_tmp_inf:
        if row[0] == "H":
            fields2 = row[9].split("_")
            if fields2[1] in sh_ucl_dict:
                if not row[3] == "100.0":
                    if sh_ucl_dict[fields2[1]] in ucl_seq_dict:
                        ucl_seq_dict[sh_ucl_dict[fields2[1]]] = ucl_seq_dict[sh_ucl_dict[fields2[1]]] + "," + row[8]
                    else:
                        ucl_seq_dict[sh_ucl_dict[fields2[1]]] = row[8]
                else:
                    o.write(f"{row[8]}\t{row[9]}\n")
                    exact_match_counter += 1
            else:
                logging.info(f"COMP\tUCL for {fields2[0]} not found.")

# print out mapping table
with open(compounds_file, "w") as o:
    for key, value in ucl_seq_dict.items():
        o.write(f"{key}\t{value}\n")
        compound_file = user_dir / "compounds" / f"{key}.fas"
        with open(compound_file, "w") as c:
            seqs = value.split(",")
            for i, val in enumerate(seqs):
                c.write(f">{val}\n")
                c.write(f"{uniques[val]}\n")

            seqs2 = seq_ucl_dict[key].split(",")
            for i, val in enumerate(seqs2):
                c.write(f">{val}\n")
                c.write(f"{original[val]}\n")
