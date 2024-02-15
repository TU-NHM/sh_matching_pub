import argparse
import csv
import os
from pathlib import Path

parser = argparse.ArgumentParser(description="Script to calculate single common taxonomy for all matches")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
matches_dir = user_dir / "matches"
infile = matches_dir / f"matches_out_all_v9.csv"
outfile = matches_dir / f"matches_out_taxonomy_v9.csv"

with open(infile) as bh, open(outfile, "w") as o:
    # print header
    o.write("seq_id" + "\t" + "seq_name" + "\t" + "common_name_selection_status" + "\t" + "common_taxonomy" + "\t" + "common_rank" + "\n")
    dataReader = csv.reader(bh, delimiter="\t")
    row_ct = 0
    for row in dataReader:
        # Indata:
        # seq_id_tmp [0]
        # seq_accno [1]
        # status (3.0) [2]
        # SH code (3.0) [3]
        # SH/compound taxonomy (3.0) [4]
        # status (2.5) [5]
        # SH code (2.5) [6]
        # SH/compound taxonomy (2.5) [7]
        # status (2.0) [8]
        # SH code (2.0) [9]
        # SH/compound taxonomy (2.0) [10]
        # status (1.5) [11]
        # SH code (1.5) [12]
        # SH/compound taxonomy (1.5) [13]
        # status (1.0) [14]
        # SH code (1.0) [15]
        # SH/compound taxonomy (1.0) [16]
        # status (0.5) [17]
        # SH code (0.5) [18]
        # SH/compound taxonomy (0.5) [19]
        # compound_cl_code (0.5) [20]
        # Compound taxonomy (0.5) [21]
        # Matched sequence [22]
        # Similarity percentage [23]
        row_ct += 1
        if row_ct > 1:
            common_taxonomy = ""
            common_rank = 0

            th_taxonomy_dict = {}
            th_taxonomy_dict[1] = {}
            th_taxonomy_dict[2] = {}
            th_taxonomy_dict[3] = {}
            th_taxonomy_dict[4] = {}
            th_taxonomy_dict[5] = {}
            th_taxonomy_dict[6] = {}

            lineage_05_arr = row[19].split(";")
            check_sp_arr_05 = None
            if len(lineage_05_arr) == 7 and not lineage_05_arr[6] == "s__unspecified":
                check_sp_arr_05 = lineage_05_arr[6][3:].split("_")
                if not check_sp_arr_05[-1] == "sp":
                    th_taxonomy_dict[6]["spe"] = lineage_05_arr[6][3:]
            if len(lineage_05_arr) > 5 and not lineage_05_arr[5] == "g__unspecified":
                check_gen_arr_05 = lineage_05_arr[5][3:].split("_")
                if len(check_gen_arr_05) == 1:
                    th_taxonomy_dict[6]["gen"] = lineage_05_arr[5][3:]
            if len(lineage_05_arr) > 4 and not lineage_05_arr[4] == "f__unspecified":
                th_taxonomy_dict[6]["fam"] = lineage_05_arr[4][3:]
            if len(lineage_05_arr) > 3 and not lineage_05_arr[3] == "o__unspecified":
                th_taxonomy_dict[6]["ord"] = lineage_05_arr[3][3:]
            if len(lineage_05_arr) > 2 and not lineage_05_arr[2] == "c__unspecified":
                th_taxonomy_dict[6]["cls"] = lineage_05_arr[2][3:]
            if len(lineage_05_arr) > 1 and not lineage_05_arr[1] == "p__unspecified":
                th_taxonomy_dict[6]["phy"] = lineage_05_arr[1][3:]
            if len(lineage_05_arr) > 0 and not lineage_05_arr[0] == "k__unspecified":
                th_taxonomy_dict[6]["kgd"] = lineage_05_arr[0][3:]
            
            lineage_10_arr = row[16].split(";")
            if len(lineage_10_arr) == 7 and not lineage_10_arr[6] == "s__unspecified":
                th_taxonomy_dict[3]["spe"] = lineage_10_arr[6][3:]
            if len(lineage_10_arr) > 5 and not lineage_10_arr[5] == "g__unspecified":
                th_taxonomy_dict[3]["gen"] = lineage_10_arr[5][3:]
            if len(lineage_10_arr) > 4 and not lineage_10_arr[4] == "f__unspecified":
                th_taxonomy_dict[3]["fam"] = lineage_10_arr[4][3:]
            if len(lineage_10_arr) > 3 and not lineage_10_arr[3] == "o__unspecified":
                th_taxonomy_dict[3]["ord"] = lineage_10_arr[3][3:]
            if len(lineage_10_arr) > 2 and not lineage_10_arr[2] == "c__unspecified":
                th_taxonomy_dict[3]["cls"] = lineage_10_arr[2][3:]
            if len(lineage_10_arr) > 1 and not lineage_10_arr[1] == "p__unspecified":
                th_taxonomy_dict[3]["phy"] = lineage_10_arr[1][3:]
            if len(lineage_10_arr) > 0 and not lineage_10_arr[0] == "k__unspecified":
                th_taxonomy_dict[3]["kgd"] = lineage_10_arr[0][3:]

            lineage_15_arr = row[13].split(";")
            if len(lineage_15_arr) == 7 and not lineage_15_arr[6] == "s__unspecified":
                th_taxonomy_dict[5]["spe"] = lineage_15_arr[6][3:]
            if len(lineage_15_arr) > 5 and not lineage_15_arr[5] == "g__unspecified":
                th_taxonomy_dict[5]["gen"] = lineage_15_arr[5][3:]
            if len(lineage_15_arr) > 4 and not lineage_15_arr[4] == "f__unspecified":
                th_taxonomy_dict[5]["fam"] = lineage_15_arr[4][3:]
            if len(lineage_15_arr) > 3 and not lineage_15_arr[3] == "o__unspecified":
                th_taxonomy_dict[5]["ord"] = lineage_15_arr[3][3:]
            if len(lineage_15_arr) > 2 and not lineage_15_arr[2] == "c__unspecified":
                th_taxonomy_dict[5]["cls"] = lineage_15_arr[2][3:]
            if len(lineage_15_arr) > 1 and not lineage_15_arr[1] == "p__unspecified":
                th_taxonomy_dict[5]["phy"] = lineage_15_arr[1][3:]
            if len(lineage_15_arr) > 0 and not lineage_15_arr[0] == "k__unspecified":
                th_taxonomy_dict[5]["kgd"] = lineage_15_arr[0][3:]

            lineage_20_arr = row[10].split(";")
            if len(lineage_20_arr) == 7 and not lineage_20_arr[6] == "s__unspecified":
                th_taxonomy_dict[2]["spe"] = lineage_20_arr[6][3:]
            if len(lineage_20_arr) > 5 and not lineage_20_arr[5] == "g__unspecified":
                th_taxonomy_dict[2]["gen"] = lineage_20_arr[5][3:]
            if len(lineage_20_arr) > 4 and not lineage_20_arr[4] == "f__unspecified":
                th_taxonomy_dict[2]["fam"] = lineage_20_arr[4][3:]
            if len(lineage_20_arr) > 3 and not lineage_20_arr[3] == "o__unspecified":
                th_taxonomy_dict[2]["ord"] = lineage_20_arr[3][3:]
            if len(lineage_20_arr) > 2 and not lineage_20_arr[2] == "c__unspecified":
                th_taxonomy_dict[2]["cls"] = lineage_20_arr[2][3:]
            if len(lineage_20_arr) > 1 and not lineage_20_arr[1] == "p__unspecified":
                th_taxonomy_dict[2]["phy"] = lineage_20_arr[1][3:]
            if len(lineage_20_arr) > 0 and not lineage_20_arr[0] == "k__unspecified":
                th_taxonomy_dict[2]["kgd"] = lineage_20_arr[0][3:]
            
            lineage_25_arr = row[7].split(";")
            if len(lineage_25_arr) == 7 and not lineage_25_arr[6] == "s__unspecified":
                th_taxonomy_dict[4]["spe"] = lineage_25_arr[6][3:]
            if len(lineage_25_arr) > 5 and not lineage_25_arr[5] == "g__unspecified":
                th_taxonomy_dict[4]["gen"] = lineage_25_arr[5][3:]
            if len(lineage_25_arr) > 4 and not lineage_25_arr[4] == "f__unspecified":
                th_taxonomy_dict[4]["fam"] = lineage_25_arr[4][3:]
            if len(lineage_25_arr) > 3 and not lineage_25_arr[3] == "o__unspecified":
                th_taxonomy_dict[4]["ord"] = lineage_25_arr[3][3:]
            if len(lineage_25_arr) > 2 and not lineage_25_arr[2] == "c__unspecified":
                th_taxonomy_dict[4]["cls"] = lineage_25_arr[2][3:]
            if len(lineage_25_arr) > 1 and not lineage_25_arr[1] == "p__unspecified":
                th_taxonomy_dict[4]["phy"] = lineage_25_arr[1][3:]
            if len(lineage_25_arr) > 0 and not lineage_25_arr[0] == "k__unspecified":
                th_taxonomy_dict[4]["kgd"] = lineage_25_arr[0][3:]
            
            lineage_30_arr = row[4].split(";")
            if len(lineage_30_arr) == 7 and not lineage_30_arr[6] == "s__unspecified":
                th_taxonomy_dict[1]["spe"] = lineage_30_arr[6][3:]
            if len(lineage_30_arr) > 5 and not lineage_30_arr[5] == "g__unspecified":
                th_taxonomy_dict[1]["gen"] = lineage_30_arr[5][3:]
            if len(lineage_30_arr) > 4 and not lineage_30_arr[4] == "f__unspecified":
                th_taxonomy_dict[1]["fam"] = lineage_30_arr[4][3:]
            if len(lineage_30_arr) > 3 and not lineage_30_arr[3] == "o__unspecified":
                th_taxonomy_dict[1]["ord"] = lineage_30_arr[3][3:]
            if len(lineage_30_arr) > 2 and not lineage_30_arr[2] == "c__unspecified":
                th_taxonomy_dict[1]["cls"] = lineage_30_arr[2][3:]
            if len(lineage_30_arr) > 1 and not lineage_30_arr[1] == "p__unspecified":
                th_taxonomy_dict[1]["phy"] = lineage_30_arr[1][3:]
            if len(lineage_30_arr) > 0 and not lineage_30_arr[0] == "k__unspecified":
                th_taxonomy_dict[1]["kgd"] = lineage_30_arr[0][3:]
            
            conflict_flag = False
            common_anc_taxonomy = ""
            if "spe" in th_taxonomy_dict[6]:
                if not check_sp_arr_05[-1] == "sp":
                    for item in th_taxonomy_dict:
                        if len(th_taxonomy_dict[item]["kgd"].split("_")) == 1:
                            if not th_taxonomy_dict[6]["kgd"] == th_taxonomy_dict[item]["kgd"]:
                                conflict_flag = True
                    if conflict_flag == True:
                        # common_anc_taxonomy = "k__conflict_on_kgd_level"
                        common_anc_taxonomy = "k__Eukaryota_kgd_Incertae_sedis"

                    if conflict_flag == False:
                        common_anc_taxonomy = th_taxonomy_dict[6]["kgd"]
                        common_rank = 1
                        for item in th_taxonomy_dict:
                            if len(th_taxonomy_dict[item]["phy"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["phy"] == th_taxonomy_dict[item]["phy"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";p__conflict_on_phy_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"]
                        common_rank = 2
                        for item in th_taxonomy_dict:
                            if len(th_taxonomy_dict[item]["cls"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["cls"] == th_taxonomy_dict[item]["cls"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";c__conflict_on_cls_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"]
                        common_rank = 3
                        for item in th_taxonomy_dict:
                            if len(th_taxonomy_dict[item]["ord"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["ord"] == th_taxonomy_dict[item]["ord"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";o__conflict_on_ord_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"]
                        common_rank = 4
                        for item in th_taxonomy_dict:
                            if len(th_taxonomy_dict[item]["fam"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["fam"] == th_taxonomy_dict[item]["fam"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";f__conflict_on_fam_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"] + ";f__" + th_taxonomy_dict[6]["fam"]
                        common_rank = 5
                        for item in th_taxonomy_dict:
                            if len(th_taxonomy_dict[item]["gen"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["gen"] == th_taxonomy_dict[item]["gen"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";g__conflict_on_gen_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"] + ";f__" + th_taxonomy_dict[6]["fam"] + ";g__" + th_taxonomy_dict[6]["gen"]
                        common_rank = 6
                        
                    if conflict_flag == False:
                        # case-1. if 0.5% species present and higher taxa not in conflict with 1.0-3.0 level taxonomy:
                        #         -> use 0.5% taxonomy:species
                        common_taxonomy = "case-1\t" + row[19] + "\t" + str(common_rank)
                    else:
                        # case-2. elif 0.5% species present but higher taxa in conflict with 1.0-3.0 level taxonomy:
                        #         -> use least common ancestor taxonomy
                        common_taxonomy = "case-2\t" + common_anc_taxonomy + "\t" + str(common_rank)
            elif "gen" in th_taxonomy_dict[6]:
                if len(check_gen_arr_05) == 1:
                    for item in th_taxonomy_dict:
                        if "kgd" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["kgd"].split("_")) == 1:
                            if not th_taxonomy_dict[6]["kgd"] == th_taxonomy_dict[item]["kgd"]:
                                conflict_flag = True
                    if conflict_flag == True:
                        # common_anc_taxonomy = "k__conflict_on_kgd_level"
                        common_anc_taxonomy = "k__Eukaryota_kgd_Incertae_sedis"

                    if conflict_flag == False:
                        common_anc_taxonomy = th_taxonomy_dict[6]["kgd"]
                        common_rank = 1
                        for item in th_taxonomy_dict:
                            if "phy" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["phy"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["phy"] == th_taxonomy_dict[item]["phy"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";p__conflict_on_phy_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"]
                        common_rank = 2
                        for item in th_taxonomy_dict:
                            if "cls" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["cls"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["cls"] == th_taxonomy_dict[item]["cls"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";c__conflict_on_cls_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"]
                        common_rank = 3
                        for item in th_taxonomy_dict:
                            if "ord" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["ord"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["ord"] == th_taxonomy_dict[item]["ord"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";o__conflict_on_ord_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"]
                        common_rank = 4
                        for item in th_taxonomy_dict:
                            if "fam" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["fam"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["fam"] == th_taxonomy_dict[item]["fam"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";f__conflict_on_fam_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"] + ";f__" + th_taxonomy_dict[6]["fam"]
                        common_rank = 5
                        for item in th_taxonomy_dict:
                            if "gen" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["gen"].split("_")) == 1:
                                if not th_taxonomy_dict[6]["gen"] == th_taxonomy_dict[item]["gen"]:
                                    conflict_flag = True
                        # if conflict_flag == True:
                        #     common_anc_taxonomy = common_anc_taxonomy + ";g__conflict_on_gen_level"

                    if conflict_flag == False:
                        common_anc_taxonomy = "k__" + th_taxonomy_dict[6]["kgd"] + ";p__" + th_taxonomy_dict[6]["phy"] + ";c__" + th_taxonomy_dict[6]["cls"] + ";o__" + th_taxonomy_dict[6]["ord"] + ";f__" + th_taxonomy_dict[6]["fam"] + ";g__" + th_taxonomy_dict[6]["gen"]
                        common_rank = 6

                    if conflict_flag == False:
                        # case-3. elif 0.5% genus present and higher taxa not in conflict with 1.0-3.0% level taxonomy:
                        #         -> use 0.5% taxonomy:genus
                        common_taxonomy = "case-3\t" + common_anc_taxonomy + "\t" + str(common_rank)
                    else:
                        # case-4. elif 0.5% genus present and genus or higher taxa in conflict with 1.0-3.0% levels:
                        #         -> use least common ancestor taxonomy
                        common_taxonomy = "case-4\t" + common_anc_taxonomy + "\t" + str(common_rank)
            else:
                common_th_taxonomy_dict = {}
                for item in th_taxonomy_dict:
                    if "kgd" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["kgd"].split("_")) == 1:
                        if "kgd" in common_th_taxonomy_dict:
                            if not common_th_taxonomy_dict["kgd"] == th_taxonomy_dict[item]["kgd"]:
                                conflict_flag = True
                        else:
                            common_th_taxonomy_dict["kgd"] = th_taxonomy_dict[item]["kgd"]
                if conflict_flag == True:
                    del common_th_taxonomy_dict["kgd"]

                if conflict_flag == False:
                    if "kgd" in common_th_taxonomy_dict:
                        common_rank = 1
                        for item in th_taxonomy_dict:
                            if "phy" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["phy"].split("_")) == 1:
                                if "phy" in common_th_taxonomy_dict:
                                    if not common_th_taxonomy_dict["phy"] == th_taxonomy_dict[item]["phy"]:
                                        conflict_flag = True
                                else:
                                    common_th_taxonomy_dict["phy"] = th_taxonomy_dict[item]["phy"]
                        if conflict_flag == True:
                            del common_th_taxonomy_dict["phy"]

                if conflict_flag == False:
                    if "phy" in common_th_taxonomy_dict:
                        common_rank = 2
                        for item in th_taxonomy_dict:
                            if "cls" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["cls"].split("_")) == 1:
                                if "cls" in common_th_taxonomy_dict:
                                    if not common_th_taxonomy_dict["cls"] == th_taxonomy_dict[item]["cls"]:
                                        conflict_flag = True
                                else:
                                    common_th_taxonomy_dict["cls"] = th_taxonomy_dict[item]["cls"]
                        if conflict_flag == True:
                            del common_th_taxonomy_dict["cls"]

                if conflict_flag == False:
                    if "cls" in common_th_taxonomy_dict:
                        common_rank = 3
                        for item in th_taxonomy_dict:
                            if "ord" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["ord"].split("_")) == 1:
                                if "ord" in common_th_taxonomy_dict:
                                    if not common_th_taxonomy_dict["ord"] == th_taxonomy_dict[item]["ord"]:
                                        conflict_flag = True
                                else:
                                    common_th_taxonomy_dict["ord"] = th_taxonomy_dict[item]["ord"]
                        if conflict_flag == True:
                            del common_th_taxonomy_dict["ord"]

                if conflict_flag == False:
                    if "ord" in common_th_taxonomy_dict:
                        common_rank = 4
                        for item in th_taxonomy_dict:
                            if "fam" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["fam"].split("_")) == 1:
                                if "fam" in common_th_taxonomy_dict:
                                    if not common_th_taxonomy_dict["fam"] == th_taxonomy_dict[item]["fam"]:
                                        conflict_flag = True
                                else:
                                    common_th_taxonomy_dict["fam"] = th_taxonomy_dict[item]["fam"]
                        if conflict_flag == True:
                            del common_th_taxonomy_dict["fam"]

                if conflict_flag == False:
                    if "fam" in common_th_taxonomy_dict:
                        common_rank = 5
                        for item in th_taxonomy_dict:
                            if "gen" in th_taxonomy_dict[item] and len(th_taxonomy_dict[item]["gen"].split("_")) == 1:
                                if "gen" in common_th_taxonomy_dict:
                                    if not common_th_taxonomy_dict["gen"] == th_taxonomy_dict[item]["gen"]:
                                        conflict_flag = True
                                else:
                                    common_th_taxonomy_dict["gen"] = th_taxonomy_dict[item]["gen"]
                        if conflict_flag == True:
                            del common_th_taxonomy_dict["gen"]

                if conflict_flag == False:
                    if "gen" in common_th_taxonomy_dict:
                        common_rank = 6

                # case-5. else:
                #     -> use least common ancestor taxonomy based on the taxonomy of 0.5-3.0% SHs
                if "kgd" in common_th_taxonomy_dict and not common_th_taxonomy_dict["kgd"] == "unidentified" and not common_th_taxonomy_dict["kgd"] == "unspecified":
                    common_anc_taxonomy = "k__" + common_th_taxonomy_dict["kgd"]
                    if "phy" in common_th_taxonomy_dict:
                        common_anc_taxonomy = common_anc_taxonomy + ";p__" + common_th_taxonomy_dict["phy"]
                        if "cls" in common_th_taxonomy_dict:
                            common_anc_taxonomy = common_anc_taxonomy + ";c__" + common_th_taxonomy_dict["cls"]
                            if "ord" in common_th_taxonomy_dict:
                                common_anc_taxonomy = common_anc_taxonomy + ";o__" + common_th_taxonomy_dict["ord"]
                                if "fam" in common_th_taxonomy_dict:
                                    common_anc_taxonomy = common_anc_taxonomy + ";f__" + common_th_taxonomy_dict["fam"]
                                    if "gen" in common_th_taxonomy_dict:
                                        common_anc_taxonomy = common_anc_taxonomy + ";g__" + common_th_taxonomy_dict["gen"]
                else:
                    common_anc_taxonomy = "k__Eukaryota_kgd_Incertae_sedis"
                common_taxonomy = "case-5\t" + common_anc_taxonomy + "\t" + str(common_rank)
            # case-6. if SH taxon name is on higher level than compound taxon name:
            #     -> use compound taxonomy instead
            compound_rank = 0
            compound_taxonomy = ""
            compound_taxonomy_arr = row[21].split(";")
            if compound_taxonomy_arr[0] and len(compound_taxonomy_arr[0][3:].split("_")) == 1 and not compound_taxonomy_arr[0] == "k__unidentified" and not compound_taxonomy_arr[0] == "k__unspecified":
                compound_taxonomy = compound_taxonomy_arr[0]
                compound_rank = 1
                if compound_taxonomy_arr[1] and len(compound_taxonomy_arr[1][3:].split("_")) == 1 and not compound_taxonomy_arr[1] == "p__unidentified" and not compound_taxonomy_arr[1] == "p__unspecified":
                    compound_taxonomy = compound_taxonomy + ";" + compound_taxonomy_arr[1]
                    compound_rank = 2
                    if compound_taxonomy_arr[2] and len(compound_taxonomy_arr[2][3:].split("_")) == 1 and not compound_taxonomy_arr[2] == "c__unidentified" and not compound_taxonomy_arr[2] == "c__unspecified":
                        compound_taxonomy = compound_taxonomy + ";" + compound_taxonomy_arr[2]
                        compound_rank = 3
                        if compound_taxonomy_arr[3] and len(compound_taxonomy_arr[3][3:].split("_")) == 1 and not compound_taxonomy_arr[3] == "o__unidentified" and not compound_taxonomy_arr[3] == "o__unspecified":
                            compound_taxonomy = compound_taxonomy + ";" + compound_taxonomy_arr[3]
                            compound_rank = 4

            if compound_rank > common_rank and conflict_flag == False:
                common_taxonomy = "case-6\t" + compound_taxonomy + "\t" + str(compound_rank)
            
            o.write(row[0] + "\t" + row[1] + "\t" + common_taxonomy + "\n")
