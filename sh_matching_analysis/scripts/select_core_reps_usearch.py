import argparse
import csv
import glob
import logging
import os
import sys

from Bio import SeqIO
from pathlib import Path

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description="Script to select 0.5 percent RepS from usearch calc_distm_+cluster_aggd clustering")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
cl_list_file = user_dir / "clusters_pre" / "tmp.txt"
singl_list_file = user_dir / "clusters_pre" / "singletons.txt"
reps_outfile = user_dir / "core_reps_pre.fasta"
infile_iupac = user_dir / "iupac_out_vsearch.fasta"
seq_mappings_file = user_dir / "seq_mappings.txt"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

# Go through usearch clustering results, print out RepS
global_counter_995 = 0
global_counter_singleton = 0

# read in all cleaned sequences
iupac_full_seqs_dict = {}
with open(infile_iupac, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        iupac_full_seqs_dict[record.id] = record.seq

# print out reps for core dataset
with open(reps_outfile, "w") as reps, open(seq_mappings_file, "a") as m:
    # Clusters
    with open(cl_list_file, "r") as f:
        dataReader = csv.reader(f, delimiter=",")
        for row in dataReader:
            # go through 995 clusters
            ucl_folder = row[0] + "_folder"
            folder_dir = user_dir / "clusters_pre" / "clusters" / ucl_folder
            if os.path.isdir(folder_dir):
                glob_match_90 = f"{folder_dir}/clusters/Cluster*"
                file_list_90 = glob.glob(glob_match_90)
                for file_90 in file_list_90:
                    ucl_code_90 = file_90.split("/")[-1]
                    tmp_us_file_90 = f"{folder_dir}/calc_distm_out/{ucl_code_90}_out_005"
                    cl_count_995_dict = {}
                    with open(tmp_us_file_90, "r") as f2:
                        dataReader2 = csv.reader(f2, delimiter="\t")
                        for row2 in dataReader2:
                            if row2[0] in cl_count_995_dict:
                                # keep discarded sequences' mappings to RepS
                                m.write(cl_count_995_dict[row2[0]] + "," + row2[1] + "\n")
                            else:
                                global_counter_995 += 1
                                cl_count_995_dict[row2[0]] = row2[1]
                                rep_seqs = row2[1]
                                reps.write(">" + rep_seqs + "\n" + str(iupac_full_seqs_dict[rep_seqs]) + "\n")

                # process singletons folder here as well
                glob_match_s_90 = f"{folder_dir}/singletons/Singleton*"
                file_list_s_90 = glob.glob(glob_match_s_90)
                for file_s_90 in file_list_s_90:
                    ucl_code_90 = file_s_90.split("/")[-1]
                    with open(file_s_90, "r") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            name = record.id
                            seq = record.seq
                            reps.write(">" + str(name) + "\n" + str(seq) + "\n")
            else:
                tmp_cl_us_out_995_name = row[0] + "_out_005"
                tmp_cl_us_out_995_file = user_dir / "clusters_pre" / "clusters" / "calc_distm_out" / tmp_cl_us_out_995_name
                cl_count_995_dict = {}
                with open(tmp_cl_us_out_995_file, "r") as f2:
                    dataReader2 = csv.reader(f2, delimiter="\t")
                    for row2 in dataReader2:
                        if row2[0] in cl_count_995_dict:
                            # keep discarded sequences' mappings to RepS
                            m.write(cl_count_995_dict[row2[0]] + "," + row2[1] + "\n")
                        else:
                            global_counter_995 += 1
                            cl_count_995_dict[row2[0]] = row2[1]
                            rep_seqs = row2[1]
                            reps.write(">" + rep_seqs + "\n" + str(iupac_full_seqs_dict[rep_seqs]) + "\n")

    # Singletons
    with open(singl_list_file, "r") as f:
        dataReader = csv.reader(f, delimiter=",")
        for row in dataReader:
            global_counter_singleton += 1
            # go through singletons
            tmp_singl_name = row[0]
            tmp_singl_file = user_dir / "clusters_pre" / "singletons" / tmp_singl_name
            with open(tmp_singl_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    name = record.id
                    seq = record.seq
                    reps.write(">" + str(name) + "\n" + str(seq) + "\n")

logging.info("No of 99.5% SHs - " + str(global_counter_995))
logging.info("No of global singletons - " + str(global_counter_singleton))
