import argparse
import csv
import os
from pathlib import Path

parser = argparse.ArgumentParser(description="Script to merge parsed matches output into one file")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
matches_dir = user_dir / "matches"
outfile = matches_dir / f"matches_out_all.csv"

## TODO: Add 995 as well once ref. data is updated
threshold_list = ["97", "975", "98", "985", "99"]
one_line_str_dict = {}

with open(outfile, "w") as o:
    for thld in threshold_list:
        matches_file = matches_dir / f"matches_out_{thld}.csv"
        with open(matches_file) as m:
            dataReader_m = csv.reader(m, delimiter="\t")
            for row in dataReader_m:
                if thld == "97":
                    one_line_str_dict[row[0]] = ""
                    ## seq_id_tmp   seq_accno   status (3.0)    SH code (3.0)   SH/compound taxonomy (3.0)
                    for i in range(5):
                        one_line_str_dict[row[0]] += row[i] + "\t"
                elif thld == "99":
                    ## status (1.0) SH code (1.0)   SH/compound taxonomy (1.0)    compound_cl_code (1.0)    Compound taxonomy (1.0)
                    one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + row[4] + "\t" + row[5] + "\t" + row[6] + "\n"
                else:
                    ## status (2.0) SH code (2.0)   SH/compound taxonomy (2.0)
                    one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + row[4] + "\t"

        matches_1_file = matches_dir / f"matches_1_out_{thld}.csv"
        with open(matches_1_file) as m1:
            dataReader_m1 = csv.reader(m1, delimiter="\t")
            row_ct = 0
            for row in dataReader_m1:
                row_ct += 1
                if row_ct > 1:
                    if thld == "97":
                        one_line_str_dict[row[0]] = ""
                        ## seq_id_tmp  seq_accno   status (3.0)    SH code (3.0)
                        for i in range(4):
                            one_line_str_dict[row[0]] += row[i] + "\t"
                        one_line_str_dict[row[0]] += "\t"
                    elif thld == "99":
                        ## status (1.0)    SH code (1.0)   compound_cl_code (1.0)
                        one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + "\t" + row[4] + "\n"
                    else:
                        ## status (3.0)    SH code (3.0)
                        one_line_str_dict[row[0]] += row[2] + "\t" + row[3] + "\t" + "\t"
    for line in one_line_str_dict:
        o.write(one_line_str_dict[line])
