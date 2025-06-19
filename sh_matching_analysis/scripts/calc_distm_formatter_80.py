import argparse
import csv
# import logging
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description="Script to run usearch single-linkage clustering for 80 percent clusters")
parser.add_argument("run_id", help="Need run id in numeric format!")
parser.add_argument("name", help="Need cluster name!")
args = parser.parse_args()

# read in args
run_id = args.run_id
name = args.name
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")

# infiles
uclust_dir = user_dir / "compounds"
out_dir = uclust_dir / "calc_distm_out"

usearch_program = "/sh_matching/programs/usearch"

code = name
code_url = uclust_dir / code
mx_code = code + "_mx_03"
mx_code_url = out_dir / mx_code
out_code_03 = code + "_out_03"
out_code_url_03 = out_dir / out_code_03
out_code_025 = code + "_out_025"
out_code_url_025 = out_dir / out_code_025
out_code_02 = code + "_out_02"
out_code_url_02 = out_dir / out_code_02
out_code_015 = code + "_out_015"
out_code_url_015 = out_dir / out_code_015
out_code_01 = code + "_out_01"
out_code_url_01 = out_dir / out_code_01
out_code_005 = code + "_out_005"
out_code_url_005 = out_dir / out_code_005

# usearch -calc_distmx ClusterX -tabbedout mx_03.txt -maxdist 0.03 -threads 8
usearch_cmd_1 = subprocess.run([usearch_program, "-calc_distmx", code_url, "-tabbedout", mx_code_url, "-maxdist", "0.03", "-threads", "8", "-gapopen", "1.0I/0.0E", "-gapext", "1.0I/0.0E"], stdout=subprocess.DEVNULL)

# usearch -cluster_aggd mx_03.txt -clusterout clusters.txt -id 0.97 -linkage min
usearch_cmd_2 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_03, "-id", "0.97", "-linkage", "min"], stdout=subprocess.DEVNULL)
usearch_cmd_3 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_025, "-id", "0.975", "-linkage", "min"], stdout=subprocess.DEVNULL)
usearch_cmd_4 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_02, "-id", "0.98", "-linkage", "min"], stdout=subprocess.DEVNULL)
usearch_cmd_5 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_015, "-id", "0.985", "-linkage", "min"], stdout=subprocess.DEVNULL)
usearch_cmd_6 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_01, "-id", "0.99", "-linkage", "min"], stdout=subprocess.DEVNULL)
usearch_cmd_7 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_005, "-id", "0.995", "-linkage", "min"], stdout=subprocess.DEVNULL)

rm_cmd_1 = subprocess.run(["rm", str(mx_code_url)], stdout=subprocess.DEVNULL)
