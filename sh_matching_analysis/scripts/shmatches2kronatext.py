import argparse
import csv
import os
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description="Script to output HTML output for matches")
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
matches_file = matches_dir / f"matches_out_{threshold}.csv"
outfile = user_dir / f"krona_{threshold}.txt"

# dicts for taxonomy div
seq_taxonomy_dict = {}
seq_taxon_count = {}

with open(outfile, "w") as o, open(matches_file) as m:
    dataReader_m = csv.reader(m, delimiter="\t")
    row_count = 0
    for row in dataReader_m:
        row_count += 1
        if row_count > 1:
            sh_taxonomy = row[4]

            # go through records' taxonomy info, populate dicts
            if sh_taxonomy and not sh_taxonomy == "":
                taxa_dict = sh_taxonomy.split(";")
                # kingdom
                if taxa_dict[0]:
                    if taxa_dict[0] in seq_taxonomy_dict:
                        seq_taxon_count[taxa_dict[0]] += 1
                    else:
                        seq_taxonomy_dict[taxa_dict[0]] = {}
                        seq_taxon_count[taxa_dict[0]] = 1
                    # phylum
                    if taxa_dict[1]:
                        tax_name_p = taxa_dict[0] + ";" + taxa_dict[1]
                        if taxa_dict[1] in seq_taxonomy_dict[taxa_dict[0]]:
                            seq_taxon_count[tax_name_p] += 1
                        else:
                            seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]] = {}
                            seq_taxon_count[tax_name_p] = 1
                        # class
                        if taxa_dict[2]:
                            tax_name_c = tax_name_p + ";" + taxa_dict[2]
                            if taxa_dict[2] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]]:
                                seq_taxon_count[tax_name_c] += 1
                            else:
                                seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]] = {}
                                seq_taxon_count[tax_name_c] = 1
                            # order
                            if taxa_dict[3]:
                                tax_name_o = tax_name_c + ";" + taxa_dict[3]
                                if taxa_dict[3] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]]:
                                    seq_taxon_count[tax_name_o] += 1
                                else:
                                    seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]] = {}
                                    seq_taxon_count[tax_name_o] = 1
                                if len(taxa_dict) > 4:
                                    # family
                                    if taxa_dict[4]:
                                        tax_name_f = tax_name_o + ";" + taxa_dict[4]
                                        if taxa_dict[4] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]:
                                            seq_taxon_count[tax_name_f] += 1
                                        else:
                                            seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]] = {}
                                            seq_taxon_count[tax_name_f] = 1
                                        # genus
                                        if taxa_dict[5]:
                                            tax_name_g = tax_name_f + ";" + taxa_dict[5]
                                            if taxa_dict[5] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]]:
                                                seq_taxon_count[tax_name_g] += 1
                                            else:
                                                seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]] = {}
                                                seq_taxon_count[tax_name_g] = 1
                                            # species
                                            if taxa_dict[6]:
                                                tax_name_s = tax_name_g + ";" + taxa_dict[6]
                                                if taxa_dict[6] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]]:
                                                    seq_taxon_count[tax_name_s] += 1
                                                else:
                                                    seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]][taxa_dict[6]] = 1
                                                    seq_taxon_count[tax_name_s] = 1
                                else:
                                    # special case for new SHs where compound taxonomy should be parsed
                                    tax_name_f = taxa_dict[0] + ";" + taxa_dict[1] + ";" + taxa_dict[2] + ";" + taxa_dict[3] + ";f__unspecified"
                                    if "f__unspecified" in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]:
                                        seq_taxon_count[tax_name_f] += 1
                                    else:
                                        seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]["f__unspecified"] = {}
                                        seq_taxon_count[tax_name_f] = 1
                                    tax_name_g = tax_name_f + ";g__unspecified"
                                    if "g__unspecified" in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]["f__unspecified"]:
                                        seq_taxon_count[tax_name_g] += 1
                                    else:
                                        seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]["f__unspecified"]["g__unspecified"] = {}
                                        seq_taxon_count[tax_name_g] = 1
                                    tax_name_s = tax_name_g + ";s__unspecified"
                                    if "s__unspecified" in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]["f__unspecified"]["g__unspecified"]:
                                        seq_taxon_count[tax_name_s] += 1
                                    else:
                                        seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]["f__unspecified"]["g__unspecified"]["s__unspecified"] = 1
                                        seq_taxon_count[tax_name_s] = 1
    # taxonomy tab
    for kgd in seq_taxonomy_dict:
        for phy in seq_taxonomy_dict[kgd]:
            for cls in seq_taxonomy_dict[kgd][phy]:
                for order in seq_taxonomy_dict[kgd][phy][cls]:
                    for fam in seq_taxonomy_dict[kgd][phy][cls][order]:
                        for gen in seq_taxonomy_dict[kgd][phy][cls][order][fam]:
                            for spec in seq_taxonomy_dict[kgd][phy][cls][order][fam][gen]:
                                tax_name_s = kgd + ";" + phy + ";" + cls + ";" + order + ";" + fam  + ";" + gen + ";" + spec
                                o.write(f"{seq_taxon_count[tax_name_s]}\t{kgd[3:]}\t{phy[3:]}\t{cls[3:]}\t{order[3:]}\t{fam[3:]}\t{gen[3:]}\t{spec[3:]}\n")
