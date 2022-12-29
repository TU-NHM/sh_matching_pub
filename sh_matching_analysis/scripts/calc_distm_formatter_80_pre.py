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
uclust_dir = user_dir / "clusters_pre" / "clusters"
out_dir = uclust_dir / "calc_distm_out"

usearch_program = "/sh_matching/programs/usearch"

code = name
code_url = uclust_dir / code
mx_code = code + "_mx_005"
mx_code_url = out_dir / mx_code
out_code_005 = code + "_out_005"
out_code_url_005 = out_dir / out_code_005

# usearch -calc_distmx ClusterX -tabbedout mx_005.txt -maxdist 0.005 -threads 8
usearch_cmd_1 = subprocess.run([usearch_program, "-calc_distmx", code_url, "-tabbedout", mx_code_url, "-maxdist", "0.005", "-threads", "8"], stdout=subprocess.DEVNULL)

# usearch -cluster_aggd mx_005.txt -clusterout clusters.txt -id 0.995 -linkage min
usearch_cmd_1 = subprocess.run([usearch_program, "-cluster_aggd", mx_code_url, "-clusterout", out_code_url_005, "-id", "0.995", "-linkage", "max"], stdout=subprocess.DEVNULL)

rm_cmd_1 = subprocess.run(["rm", str(mx_code_url)], stdout=subprocess.DEVNULL)
