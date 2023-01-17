import argparse
import csv
import os
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
matches_1_file = matches_dir / f"matches_1_out_{threshold}.csv"
outfile = matches_dir / f"matches_out_{threshold}.html"

style_tag = "<style>\nbody {background-color: #d0dadd; padding: 5px; }\n\
span, p {color: #333333; font-family: \"Helvetica Neue\", Helvetica, Roboto, Arial, sans-serif; font-style: normal; font-weight: normal; }\n\
table, th, td {border-collapse: collapse; border: 1px solid #222222;}\
th, td {padding: 2px 4px; }\
tr:nth-child(even) {background-color: #eaeff0;}\
th {text-align: left; font-size: 0.9rem; color: #333333; font-family: \"Helvetica Neue\", Helvetica, Roboto, Arial, sans-serif; font-style: normal;}\n\
td {font-size: 0.9rem; color: #333333; font-family: \"Helvetica Neue\", Helvetica, Roboto, Arial, sans-serif; font-style: normal; font-weight: normal;} \
h5 {font-size: 1.125rem; line-height: 1.4; color: #222222; margin-bottom: 0.5rem; margin-top: 0.2rem; font-family: \"Helvetica Neue\", Helvetica, Roboto, Arial, sans-serif; font-style: normal; font-weight: normal;}\n\
h4 {font-size: 1rem; color: #333333; font-family: \"Helvetica Neue\", Helvetica, Roboto, Arial, sans-serif; font-style: normal; font-weight: normal;}\n\
a {color: #337ab7}\n\
.small_text { font-size: 0.9rem; }\
</style>"

threshold_dict = {"03": "3.0", "025": "2.5", "02": "2.0", "015": "1.5", "01": "1.0", "005": "0.5"}

# dicts for statistics div
# seqs_status_hash = {}
shs_hash = {}
shs_hash["present_in"] = {}
shs_hash["new_sh_in"] = {}
shs_hash["new_singleton_in"] = {}
ex_shs_count = 0
new_shs_count = 0
new_singleton_shs_count = 0
ex_seqs_count = 0
new_seqs_count = 0
new_singleton_seqs_count = 0
ucl_hash = {}
ucl_hash["present_in"] = {}
ucl_hash["new_sh_in"] = {}
ucl_hash["new_singleton_in"] = {}
ucl_hash["total"] = {}
ex_ucl_count = 0
new_ucl_count = 0
new_singleton_ucl_count = 0
total_ucl_count = 0

# new compounds
nc_shs_hash = {}
nc_shs_hash["new_sh_in"] = {}
nc_shs_hash["new_singleton_in"] = {}
nc_ucl_hash = {}
nc_ucl_hash["new_sh_in"] = {}
nc_ucl_hash["new_singleton_in"] = {}
nc_shs_count = 0
nc_singleton_shs_count = 0
nc_seqs_count = 0
nc_singleton_seqs_count = 0
nc_ucl_count = 0
nc_singleton_ucl_count = 0

# dicts for taxonomy div
seq_taxonomy_dict = {}
seq_taxon_count = {}
seq_taxon_sh_count = {}

with open(outfile, "w") as o, open(matches_file) as m:
    # print HTML file header
    o.write(f"<!DOCTYPE html><html>\n{style_tag}\n<head>\n<title>SH matching analysis: source_{run_id}</title>\n<head>\n<body>\n")
    o.write(f"<h5>SH matching analysis: source_{run_id}</h5>\n")
    o.write(f"<span>Distance threshold: {threshold_dict[threshold]}%</span>")
    dataReader_m = csv.reader(m, delimiter="\t")
    row_count = 0
    record_count = 0
    for row in dataReader_m:
        row_count += 1
        if row_count > 1:
            record_count += 1

            status = row[2]
            sh_code = row[3]
            sh_taxonomy = row[4]
            compound_code = row[5]

            # count existing SHs
            if status == "present_in":
                if sh_code in shs_hash["present_in"]:
                    shs_hash["present_in"][sh_code] += 1
                else:
                    shs_hash["present_in"][sh_code] = 1
                    ex_shs_count += 1
                ex_seqs_count += 1
            elif status == "new_sh_in":
                if sh_code in shs_hash["new_sh_in"]:
                    shs_hash["new_sh_in"][sh_code] += 1
                else:
                    shs_hash["new_sh_in"][sh_code] = 1
                    new_shs_count += 1
                new_seqs_count += 1
            elif status == "new_singleton_in":
                if sh_code in shs_hash["new_singleton_in"]:
                    shs_hash["new_singleton_in"][sh_code] += 1
                else:
                    shs_hash["new_singleton_in"][sh_code] = 1
                    new_singleton_shs_count += 1
                new_singleton_seqs_count += 1
            # count existing compounds
            if status == "present_in":
                if compound_code in ucl_hash["present_in"]:
                    ucl_hash["present_in"][compound_code] += 1
                else:
                    ucl_hash["present_in"][compound_code] = 1
                    ex_ucl_count += 1
            elif status == "new_sh_in":
                if compound_code in ucl_hash["new_sh_in"]:
                    ucl_hash["new_sh_in"][compound_code] += 1
                else:
                    ucl_hash["new_sh_in"][compound_code] = 1
                    new_ucl_count += 1
            elif status == "new_singleton_in":
                if compound_code in ucl_hash["new_singleton_in"]:
                    ucl_hash["new_singleton_in"][compound_code] += 1
                else:
                    ucl_hash["new_singleton_in"][compound_code] = 1
                    new_singleton_ucl_count += 1
            if not compound_code in ucl_hash["total"]:
                ucl_hash["total"][compound_code] = 1
                total_ucl_count += 1
            else:
                ucl_hash["total"][compound_code] += 1

            # go through records' taxonomy info, populate dicts
            if sh_taxonomy and not sh_taxonomy == "":
                taxa_dict = sh_taxonomy.split(";")
                # kingdom
                if taxa_dict[0]:
                    if taxa_dict[0] in seq_taxonomy_dict:
                        seq_taxon_count[taxa_dict[0]] += 1
                        if sh_code in seq_taxon_sh_count[taxa_dict[0]]:
                            seq_taxon_sh_count[taxa_dict[0]][sh_code] += 1
                        else:
                            seq_taxon_sh_count[taxa_dict[0]][sh_code] = 1
                    else:
                        seq_taxonomy_dict[taxa_dict[0]] = {}
                        seq_taxon_count[taxa_dict[0]] = 1
                        seq_taxon_sh_count[taxa_dict[0]] = {}
                        seq_taxon_sh_count[taxa_dict[0]][sh_code] = 1
                    # phylum
                    if taxa_dict[1]:
                        tax_name_p = taxa_dict[0] + ";" + taxa_dict[1]
                        if taxa_dict[1] in seq_taxonomy_dict[taxa_dict[0]]:
                            seq_taxon_count[tax_name_p] += 1
                            if sh_code in seq_taxon_sh_count[tax_name_p]:
                                seq_taxon_sh_count[tax_name_p][sh_code] += 1
                            else:
                                seq_taxon_sh_count[tax_name_p][sh_code] = 1
                        else:
                            seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]] = {}
                            seq_taxon_count[tax_name_p] = 1
                            seq_taxon_sh_count[tax_name_p] = {}
                            seq_taxon_sh_count[tax_name_p][sh_code] = 1
                        # class
                        if taxa_dict[2]:
                            tax_name_c = tax_name_p + ";" + taxa_dict[2]
                            if taxa_dict[2] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]]:
                                seq_taxon_count[tax_name_c] += 1
                                if sh_code in seq_taxon_sh_count[tax_name_c]:
                                    seq_taxon_sh_count[tax_name_c][sh_code] += 1
                                else:
                                    seq_taxon_sh_count[tax_name_c][sh_code] = 1
                            else:
                                seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]] = {}
                                seq_taxon_count[tax_name_c] = 1
                                seq_taxon_sh_count[tax_name_c] = {}
                                seq_taxon_sh_count[tax_name_c][sh_code] = 1
                            # order
                            if taxa_dict[3]:
                                tax_name_o = tax_name_c + ";" + taxa_dict[3]
                                if taxa_dict[3] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]]:
                                    seq_taxon_count[tax_name_o] += 1
                                    if sh_code in seq_taxon_sh_count[tax_name_o]:
                                        seq_taxon_sh_count[tax_name_o][sh_code] += 1
                                    else:
                                        seq_taxon_sh_count[tax_name_o][sh_code] = 1
                                else:
                                    seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]] = {}
                                    seq_taxon_count[tax_name_o] = 1
                                    seq_taxon_sh_count[tax_name_o] = {}
                                    seq_taxon_sh_count[tax_name_o][sh_code] = 1
                                if len(taxa_dict) > 4:
                                    # family
                                    if taxa_dict[4]:
                                        tax_name_f = tax_name_o + ";" + taxa_dict[4]
                                        if taxa_dict[4] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]]:
                                            seq_taxon_count[tax_name_f] += 1
                                            if sh_code in seq_taxon_sh_count[tax_name_f]:
                                                seq_taxon_sh_count[tax_name_f][sh_code] += 1
                                            else:
                                                seq_taxon_sh_count[tax_name_f][sh_code] = 1
                                        else:
                                            seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]] = {}
                                            seq_taxon_count[tax_name_f] = 1
                                            seq_taxon_sh_count[tax_name_f] = {}
                                            seq_taxon_sh_count[tax_name_f][sh_code] = 1
                                        # genus
                                        if taxa_dict[5]:
                                            tax_name_g = tax_name_f + ";" + taxa_dict[5]
                                            if taxa_dict[5] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]]:
                                                seq_taxon_count[tax_name_g] += 1
                                                if sh_code in seq_taxon_sh_count[tax_name_g]:
                                                    seq_taxon_sh_count[tax_name_g][sh_code] += 1
                                                else:
                                                    seq_taxon_sh_count[tax_name_g][sh_code] = 1
                                            else:
                                                seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]] = {}
                                                seq_taxon_count[tax_name_g] = 1
                                                seq_taxon_sh_count[tax_name_g] = {}
                                                seq_taxon_sh_count[tax_name_g][sh_code] = 1
                                            # species
                                            if taxa_dict[6]:
                                                tax_name_s = tax_name_g + ";" + taxa_dict[6]
                                                if taxa_dict[6] in seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]]:
                                                    seq_taxon_count[tax_name_s] += 1
                                                    if sh_code in seq_taxon_sh_count[tax_name_s]:
                                                        seq_taxon_sh_count[tax_name_s][sh_code] += 1
                                                    else:
                                                        seq_taxon_sh_count[tax_name_s][sh_code] = 1
                                                else:
                                                    seq_taxonomy_dict[taxa_dict[0]][taxa_dict[1]][taxa_dict[2]][taxa_dict[3]][taxa_dict[4]][taxa_dict[5]][taxa_dict[6]] = 1
                                                    seq_taxon_count[tax_name_s] = 1
                                                    seq_taxon_sh_count[tax_name_s] = {}
                                                    seq_taxon_sh_count[tax_name_s][sh_code] = 1

    if matches_1_file.exists():
        with open(matches_1_file) as m1:
            # count sequences and SHs in new compound clusters
            dataReader_m1 = csv.reader(m1, delimiter="\t")
            row_count = 0
            for row in dataReader_m1:
                row_count += 1
                if row_count > 1:
                    record_count += 1

                    status = row[2]
                    sh_code = row[3]
                    compound_code = row[4]

                    # count new SHs
                    if status == "new_sh_in":
                        if sh_code in nc_shs_hash["new_sh_in"]:
                            nc_shs_hash["new_sh_in"][sh_code] += 1
                        else:
                            nc_shs_hash["new_sh_in"][sh_code] = 1
                            nc_shs_count += 1
                        nc_seqs_count += 1
                    elif status == "new_singleton_in":
                        if sh_code in nc_shs_hash["new_singleton_in"]:
                            nc_shs_hash["new_singleton_in"][sh_code] += 1
                        else:
                            nc_shs_hash["new_singleton_in"][sh_code] = 1
                            nc_singleton_shs_count += 1
                        nc_singleton_seqs_count += 1
                    # count new compounds
                    if status == "new_sh_in":
                        if compound_code in nc_ucl_hash["new_sh_in"]:
                            nc_ucl_hash["new_sh_in"][compound_code] += 1
                        else:
                            nc_ucl_hash["new_sh_in"][compound_code] = 1
                            nc_ucl_count += 1
                    elif status == "new_singleton_in":
                        if compound_code in nc_ucl_hash["new_singleton_in"]:
                            nc_ucl_hash["new_singleton_in"][compound_code] += 1
                        else:
                            nc_ucl_hash["new_singleton_in"][compound_code] = 1
                            nc_singleton_ucl_count += 1
                    if not compound_code in ucl_hash["total"]:
                        ucl_hash["total"][compound_code] = 1
                        total_ucl_count += 1
                    else:
                        ucl_hash["total"][compound_code] += 1

    # sequence statistics tab
    o.write("<h4><b>Statistics</b></h4>\n")
    o.write(f"<table>\n<tr><th>Status</th><th>#sequences</th><th>#SHs ({threshold_dict[threshold]}%)</th><th>#compounds</th></tr>\n")
    # present_in
    o.write("<tr><td>Present in existing SHs (&quot;present_in&quot;)</td><td>")
    o.write(str(ex_seqs_count))
    o.write("</td><td>")
    o.write(str(ex_shs_count))
    o.write("</td><td>")
    o.write(str(ex_ucl_count))
    o.write("</td></tr>\n")
    # new_sh_in
    o.write("<tr><td>Form new SHs in existing compound clusters (&quot;new_sh_in&quot;)</td><td>")
    o.write(str(new_seqs_count))
    o.write("</td><td>")
    o.write(str(new_shs_count))
    o.write("</td><td>")
    o.write(str(new_ucl_count))
    o.write("</td></tr>\n")
    # new_singleton_in
    o.write("<tr><td>Form new singleton SHs in existing compound clusters (&quot;new_singleton_in&quot;)</td><td>")
    o.write(str(new_singleton_seqs_count))
    o.write("</td><td>")
    o.write(str(new_singleton_shs_count))
    o.write("</td><td>")
    o.write(str(new_singleton_ucl_count))
    o.write("</td></tr>\n")
    # new compound: new_sh_in
    o.write("<tr><td>Form new SHs in new compound clusters (&quot;new_sh_in&quot;)</td><td>")
    o.write(str(nc_seqs_count))
    o.write("</td><td>")
    o.write(str(nc_shs_count))
    o.write("</td><td>")
    o.write(str(nc_ucl_count))
    o.write("</td></tr>\n")
    # new compound: new_singleton_in
    o.write("<tr><td>Form new singleton SHs in new compound clusters (&quot;new_singleton_in&quot;)</td><td>")
    o.write(str(nc_singleton_seqs_count))
    o.write("</td><td>")
    o.write(str(nc_singleton_shs_count))
    o.write("</td><td>")
    o.write(str(nc_singleton_ucl_count))
    o.write("</td></tr>\n")
    # total
    total_new_seqs_ct = ex_seqs_count + new_seqs_count + new_singleton_seqs_count + nc_seqs_count + nc_singleton_seqs_count
    total_new_shs_ct = ex_shs_count + new_shs_count + new_singleton_shs_count + nc_shs_count + nc_singleton_shs_count
    o.write("<tr><td><strong>Total (unique)</strong></td><td><strong>")
    o.write(str(total_new_seqs_ct))
    o.write("</strong></td><td><strong>")
    o.write(str(total_new_shs_ct))
    o.write("</strong></td><td><strong>")
    o.write(str(total_ucl_count))
    o.write("</strong></td></tr>\n")
    o.write("</table>\n")

    # taxonomy tab
    o.write("<br />\n")
    o.write("<h4><b>Taxonomy</b></h4>\n")
    o.write("<p class='small_text'>NB! Only DNA sequences placed into existing compound clusters are included in the following taxonomy distribution table.</p>\n")
    o.write(f"<table>\n<tr><th>Taxon name</th><th>#sequences</th><th><nobr>#SHs ({threshold_dict[threshold]}%)</nobr></th><th>List of {threshold_dict[threshold]}% SHs (sequence count). Only SHs with matching record status set to &quot;present_in&quot; are listed here.</th></tr>\n")
    # grep 'present_in' userdir/3/matches/matches_out_97.csv | awk -F '\t' '{print $6}'| awk -F ';' '{print $1}' | sort| uniq -c | sort -k1 -rn
    for kgd in seq_taxonomy_dict:
        o.write(f"<tr><td><b>{kgd}</b></td><td>{seq_taxon_count[kgd]}</td><td>{len(seq_taxon_sh_count[kgd])}</td><td></td></tr>")
        for phy in seq_taxonomy_dict[kgd]:
            tax_name_p = kgd + ";" + phy
            o.write(f"<tr><td>&emsp;{phy}</td><td>{seq_taxon_count[tax_name_p]}</td><td>{len(seq_taxon_sh_count[tax_name_p])}</td><td></td></tr>")
            for clss in seq_taxonomy_dict[kgd][phy]:
                tax_name_c = tax_name_p + ";" + clss
                o.write(f"<tr><td>&emsp;&emsp;{clss}</td><td>{seq_taxon_count[tax_name_c]}</td><td>{len(seq_taxon_sh_count[tax_name_c])}</td><td></td></tr>")
                for order in seq_taxonomy_dict[kgd][phy][clss]:
                    tax_name_o = tax_name_c + ";" + order
                    o.write(f"<tr><td>&emsp;&emsp; &emsp;{order}</td><td>{seq_taxon_count[tax_name_o]}</td><td>{len(seq_taxon_sh_count[tax_name_o])}</td><td>")
                    # if len(seq_taxonomy_dict[kgd][phy][clss][order]) == 0:
                    #     for sh in seq_taxon_sh_count[tax_name_o]:
                    #         o.write(f"<a href='https://unite.ut.ee/sh/{sh}' target='_blank'>{sh}</a> ({seq_taxon_sh_count[tax_name_o][sh]}); ")
                    o.write("</td></tr>")
                    for fam in seq_taxonomy_dict[kgd][phy][clss][order]:
                        tax_name_f = tax_name_o + ";" + fam
                        o.write(f"<tr><td>&emsp;&emsp; &emsp;&emsp;{fam}</td><td>{seq_taxon_count[tax_name_f]}</td><td>{len(seq_taxon_sh_count[tax_name_f])}</td><td></td></tr>")
                        for gen in seq_taxonomy_dict[kgd][phy][clss][order][fam]:
                            tax_name_g = tax_name_f + ";" + gen
                            o.write(f"<tr><td>&emsp;&emsp; &emsp;&emsp; &emsp;{gen}</td><td>{seq_taxon_count[tax_name_g]}</td><td>{len(seq_taxon_sh_count[tax_name_g])}</td><td></td></tr>")
                            for spe in seq_taxonomy_dict[kgd][phy][clss][order][fam][gen]:
                                tax_name_s = tax_name_g + ";" + spe
                                o.write(f"<tr><td>&emsp;&emsp; &emsp;&emsp; &emsp;&emsp;{spe}</td><td>{seq_taxon_count[tax_name_s]}</td><td>{len(seq_taxon_sh_count[tax_name_s])}</td><td>")
                                for sh in seq_taxon_sh_count[tax_name_s]:
                                    o.write(f"<a href='https://unite.ut.ee/sh/{sh}' target='_blank'>{sh}</a> ({seq_taxon_sh_count[tax_name_s][sh]}); ")
                                o.write("</td></tr>")
    o.write("</table>\n")

    o.write(f"</body>\n<html>")
