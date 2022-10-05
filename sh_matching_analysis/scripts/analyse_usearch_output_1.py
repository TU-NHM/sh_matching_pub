import argparse
import csv
import glob
import os
import re
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
matches_dir = user_dir / "matches"
matches_file = matches_dir / f"matches_1_{threshold}.txt"

prev_dict = {}

# read in results from prev version (if singleton, no need to check it again)
prev_file = None
if threshold == "025":
    prev_file = matches_dir / "matches_1_03.txt"
elif threshold == "02":
    prev_file = matches_dir / "matches_1_025.txt"
elif threshold == "015":
    prev_file = matches_dir / "matches_1_02.txt"
elif threshold == "01":
    prev_file = matches_dir / "matches_1_015.txt"
elif threshold == "005":
    prev_file = matches_dir / "matches_1_01.txt"

if prev_file:
    with open(prev_file) as p:
        dataReader = csv.reader(p, delimiter="\t")
        for row in dataReader:
            if row[2] == "singleton":
                prev_dict[row[0]] = 1

# open matches file
with open(matches_file, "w") as o:
    # go through usearch output
    print("Processing SHs at " + str(threshold) + "...")
    glob_match = f"{user_dir}/clusters/clusters/Cluster*"
    file_list = glob.glob(glob_match)
    for file in file_list:
        if re.search(r'Cluster\d{1,}$', file) is not None:
            ucl_code = file.split("/")[-1]
            # check if 90% cluster folder exists
            ucl_folder = ucl_code + "_folder"
            folder_dir = user_dir / "clusters" / "clusters" / ucl_folder
            if os.path.isdir(folder_dir):
                glob_match_90 = f"{folder_dir}/Cluster*"
                file_list_90 = glob.glob(glob_match_90)
                for file_90 in file_list_90:
                    ucl_code_90 = file_90.split("/")[-1]
                    tmp_bc_file_90 = f"{folder_dir}/calc_distm_out/{ucl_code_90}_out_{threshold}"
                    cl_contents_dict_90 = {}
                    with open(tmp_bc_file_90) as f:
                        dataReader_90 = csv.reader(f_90, delimiter="\t")
                        for row_90 in dataReader_90:
                            if row_90[0] in cl_contents_dict_90:
                                cl_contents_dict_90[row_90[0]] = cl_contents_dict_90[row_90[0]] + " " + row_90[1]
                            else:
                                cl_contents_dict_90[row_90[0]] = row_90[1]
                        for key_90 in cl_contents_dict_90:
                            row2_90 = cl_contents_dict_90[key_90].split(" ")
                            if row2_90[0] not in prev_dict:
                                if len(row2_90) == 1:
                                    o.write(f"{row2_90[0]}\t\tsingleton\t{ucl_code}\t\n")
                                else:
                                    o.write(f"{row2_90[0]}\t\tnew cluster\t{ucl_code}\t{cl_contents_dict_90[key_90]}\n")
                            else:
                                o.write(f"{row2_90[0]}\t\tsingleton\t{ucl_code}\t\n")
                # process singletons folder here as well
                glob_match_s_90 = f"{folder_dir}/clusters/singletons/Singleton*"
                file_list_s_90 = glob.glob(glob_match_s_90)
                for file_s_90 in file_list_s_90:
                    ucl_code_90 = file_s_90.split("/")[-1]
                    with open(file_s_90, "r") as handle:
                        for record in SeqIO.parse(handle, "fasta"):
                            if record.id not in prev_dict:
                                o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")
                            else:
                                o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")

            else:
                ucl_code = file.split("/")[-1]
                tmp_bc_file = f"{user_dir}/clusters/clusters/calc_distm_out/{ucl_code}_out_{threshold}"
                cl_contents_dict_80 = {}
                with open(tmp_bc_file) as f:
                    dataReader = csv.reader(f, delimiter="\t")
                    for row_80 in dataReader:
                        if row_80[0] in cl_contents_dict_80:
                            cl_contents_dict_80[row_80[0]] = cl_contents_dict_80[row_80[0]] + " " + row_80[1]
                        else:
                            cl_contents_dict_80[row_80[0]] = row_80[1]
                    for key_80 in cl_contents_dict_80:
                        row2_80 = cl_contents_dict_80[key_80].split(" ")
                        if row2_80[0] not in prev_dict:
                            if len(row2_80) == 1:
                                o.write(f"{row2_80[0]}\t\tsingleton\t{ucl_code}\t\n")
                            else:
                                joined_row = " ".join(row2_80)
                                for seq in row2_80:
                                    if not seq == "":
                                        o.write(f"{seq}\t\tnew cluster\t{ucl_code}\t{joined_row}\n")
                        else:
                            o.write(f"{row2_80[0]}\t\tsingleton\t{ucl_code}\t\n")
    print("Finished with SHs.")

    # process singletons folder
    print("Processing global singletons at " + str(threshold) + "...")
    glob_match_s = f"{user_dir}/clusters/singletons/Singleton*"
    file_list_s = glob.glob(glob_match_s)
    for file_s in file_list_s:
        ucl_code = file_s.split("/")[-1]
        with open(file_s, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id not in prev_dict:
                    o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")
                else:
                    o.write(f"{record.id}\t\tsingleton\t{ucl_code}\t\n")
    print("Finished with global singletons.\n")
