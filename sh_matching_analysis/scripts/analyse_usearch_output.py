import argparse
import csv
import glob
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description="Script to parse usearch cl_aggd output")
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
sh2compound_file = "/sh_matching/data/sh2compound_mapping.txt"
matches_dir = user_dir / "matches"
matches_file = matches_dir / f"matches_{threshold}.txt"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

# read in results from prev version (if singleton, no need to check it again)
prev_file = None
if threshold == "025":
    prev_file = matches_dir / "matches_03.txt"
elif threshold == "02":
    prev_file = matches_dir / "matches_025.txt"
elif threshold == "015":
    prev_file = matches_dir / "matches_02.txt"
elif threshold == "01":
    prev_file = matches_dir / "matches_015.txt"
elif threshold == "005":
    prev_file = matches_dir / "matches_01.txt"

prev_dict = {}
if prev_file:
    with open(prev_file) as p:
        dataReader = csv.reader(p, delimiter="\t")
        for row in dataReader:
            if row[2] == "singleton":
                prev_dict[row[0]] = 1

# get compound and SH mappings
sh_ucl_dict = {}
with open(sh2compound_file) as f:
    dataReader = csv.reader(f, delimiter="\t")
    for row in dataReader:
        sh_ucl_dict[row[0]] = row[1]

# read seq ids and their UCL belongings into hash
seq_ucl_mapping_dict = {}
glob_match = f"{user_dir}/compounds/calc_distm_out/*.fas_out_{threshold}"
file_list = glob.glob(glob_match)
for file in file_list:
    with open(file) as f:
        dataReader = csv.reader(f, delimiter="\t")
        for row in dataReader:
            seq_ucl_mapping_dict[row[1]] = file

glob_match_folders = f"{user_dir}/compounds/*_folder"
folder_list = glob.glob(glob_match_folders)
for folder in folder_list:
    compound_name = folder[:-7].replace("compounds", "compounds/calc_distm_out")
    # read in clusters
    glob_match_files = f"{folder}/calc_distm_out/*_out_{threshold}"
    file_list = glob.glob(glob_match_files)
    for file in file_list:
        with open(file) as f:
            dataReader = csv.reader(f, delimiter="\t")
            for row in dataReader:
                seq_ucl_mapping_dict[row[1]] = f"{compound_name}_out_{threshold}"
    # read in singletons
    s_list_file = f"{folder}/singletons.txt"
    with open (s_list_file, "r") as l:
        dataReader = csv.reader(l, delimiter="\t")
        for row in dataReader:
            s_file = f"{folder}/singletons/{row[0]}"
            with open(s_file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    seq_ucl_mapping_dict[str(record.id)] = f"{compound_name}_out_{threshold}"

# go through usearch output, write out BC style clustering results
tmp_file = user_dir / "compounds" / "tmp.txt"
with open(tmp_file) as t:
    dataReader = csv.reader(t, delimiter="\t")
    for row in dataReader:
        ucl_code = row[0]
        folder_dir = user_dir / "compounds" / f"{ucl_code}_folder"
        # print out clustering results in BC style
        tmp_bc_file_new = user_dir / "compounds" / "calc_distm_out" / f"{ucl_code}_out_{threshold}_bc"
        with open(tmp_bc_file_new, "w") as b_90:
            if os.path.isdir(folder_dir):
                tmp_infile_90 = folder_dir / "tmp.txt"
                with open(tmp_infile_90) as i_90:
                    dataReader_i_90 = csv.reader(i_90, delimiter="\t")
                    for row_i_90 in dataReader_i_90:
                        ucl_code_90 = row_i_90[0]
                        # read in clustering results in usearch style
                        cl_contents_dict_90 = {}
                        tmp_bc_file_90 = folder_dir / "calc_distm_out" / f"{ucl_code_90}_out_{threshold}"
                        with open(tmp_bc_file_90) as bc_90:
                            dataReader_bc_90 = csv.reader(bc_90, delimiter="\t")
                            for row_bc_90 in dataReader_bc_90:
                                if row_bc_90[0] in cl_contents_dict_90:
                                    cl_contents_dict_90[row_bc_90[0]] = cl_contents_dict_90[row_bc_90[0]] + " " + row_bc_90[1]
                                else:
                                    cl_contents_dict_90[row_bc_90[0]] = row_bc_90[1]
                                
                        for key_90 in cl_contents_dict_90:
                            b_90.write(cl_contents_dict_90[key_90] + "\n")

                # include 90% singletons here as well
                singl_infile_90 = folder_dir / "singletons.txt";
                with open(singl_infile_90) as i_90_s:
                    dataReader_i_90_s = csv.reader(i_90_s, delimiter="\t")
                    for row_i_90_s in dataReader_i_90_s:
                        file_s_90 = folder_dir / "singletons" / f"{row_i_90_s[0]}"
                        with open(file_s_90, "r") as handle:
                            for record in SeqIO.parse(handle, "fasta"):
                                b_90.write(str(record.id) + "\n")
            else:
                # normal cluster
                # read in clustering results in usearch style
                tmp_bc_file = user_dir / "compounds" / "calc_distm_out" / f"{ucl_code}_out_{threshold}"
                # print out clustering results in BC style
                cl_contents_dict_80 = {}
                with open(tmp_bc_file) as bc_in:
                    dataReader_i_80 = csv.reader(bc_in, delimiter="\t")
                    for row_i_80 in dataReader_i_80:
                        if row_i_80[0] in cl_contents_dict_80:
                            cl_contents_dict_80[row_i_80[0]] = cl_contents_dict_80[row_i_80[0]] + " " + row_i_80[1]
                        else:
                            cl_contents_dict_80[row_i_80[0]] = row_i_80[1];
                for key_80 in cl_contents_dict_80:
                    b_90.write(cl_contents_dict_80[key_80] + "\n")

# read in best matches (seq and SH)
tmp_uc_infile = user_dir / "closedref.80-best-hits.map.uc"
with open(tmp_uc_infile) as t, open(matches_file, "w") as o:
    dataReader = csv.reader(t, delimiter="\t")
    for row in dataReader:
        if row[0] == "H":
            fields = row[9].split("_", 2)
            query = row[8]
            subject = fields[0]
            subject_sh = fields[1]
            check_output = None
            if query not in prev_dict:
                o.write(query + "\t" + subject + "\t")
                tmp_folder = user_dir / "compounds" / "calc_distm_out" / f"{sh_ucl_dict[subject_sh]}.fas_out_{threshold}_bc"
                check = subprocess.run(f"grep {query} {tmp_folder} | grep {subject}", shell=True, universal_newlines=True, stdout=subprocess.PIPE)
                check_output = check.stdout.rstrip()

                if check_output and not check_output == "":
                    o.write("present" + "\t" + sh_ucl_dict[subject_sh] + "\t" + "\n")
                else:
                    bc_tmp_folder = seq_ucl_mapping_dict[query] + "_bc"
                    check = subprocess.run(f"grep {query} {bc_tmp_folder}", shell=True, universal_newlines=True, stdout=subprocess.PIPE)
                    check_output = check.stdout.rstrip()
                    if check_output and not check_output == "":
                        if check_output == query:
                            o.write("singleton\t" + sh_ucl_dict[subject_sh] + "\t" + "\n")
                        else:
                            # check if cluster contents includes only query sequences
                            check_output_list = check_output.split(" ")
                            new_cl_status = True
                            seq_match = None
                            for list_item in check_output_list:
                                if not "_" in list_item:
                                    new_cl_status = False
                                    seq_match = list_item
                                    break
                            if new_cl_status == True:
                                o.write("new cluster\t" + sh_ucl_dict[subject_sh] + "\t" + check_output + "\n")
                            else:
                                o.write("present\t" + sh_ucl_dict[subject_sh] + "\t" + "\n")
                                logging.info(f"US\tcorrect RepS not matched for {query}")
                    else:
                        o.write("missing\t\t\n")
            else:
                o.write(query + "\t" + subject + "\t")
                o.write("singleton\t" + sh_ucl_dict[subject_sh] + "\t" + "\n")
        else:
            logging.info(f"US\tno hit for {row[8]}\n")
