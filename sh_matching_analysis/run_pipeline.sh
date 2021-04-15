#!/bin/bash

# running the script - ./run_pipeline.sh <run_id>

# TODO LIST:
# 1. print all excluded sequences in a separate file together with the reason they were excluded - Done
# 2. generate output files

if [ -z "$1" ]
    then
        echo "No run_id argument supplied!"
        exit
fi

# get run id
run_id=$1
region=$2

if [ "$region" != "its2" ] && [ "$region" != "itsfull" ]; then
  echo "Setting region to itsfull"
  region="itsfull"
fi

# get working directory
pwd=$(pwd)

script_dir="/sh_matching/scripts"
data_dir="/sh_matching/data"
outdata_dir="$pwd/outdata"
user_dir="$pwd/userdir/"$run_id
program_dir="/sh_matching/programs"
infile="source_$run_id"
infile_new=$infile"_fasta"
infile_new_w_dir="$user_dir/$infile_new"

echo "Start"

echo "Starting new analysis in $user_dir"

# rm (if exists) user working dir and create as new
if [ -d "$user_dir" ]
  then
      rm -fr "$user_dir"
fi
mkdir "$user_dir"

err_log="$user_dir/err_$run_id.log"
touch "$err_log"
ex_log="$user_dir/excluded_$run_id.txt"
touch "$ex_log"

# copy indata to user's workdir
cp "$pwd/indata/$infile" "$user_dir"

# replace sequence identifiers with unique codes for the analysis
python3 "$script_dir"/replace_seq_names_w_codes.py "$run_id"

# remove duplicate sequences from user’s dataset
pushd "$user_dir"
mothur "#unique.seqs(fasta=$infile_new_w_dir)"
mv "$infile_new""unique" "$infile_new""unique_mothur"
popd

python3 "$script_dir"/reformat_mothur_output.py "$run_id"

# Extract ITS regions, first for fungi, and then all other groups.
pushd "$user_dir"
mkdir ITSx
perl "$program_dir/ITSx/ITSx" -t F -i "$infile_new""unique" -o itsx_sh_out --cpu 8 --save_regions all --partial 50 --detailed_results T -concat T -preserve T -N 1 -E 0.01 --search_eval 0.01
mv itsx_sh_out* ITSx/

mkdir ITSx_o
perl "$program_dir/ITSx/ITSx" -t A,B,C,D,E,G,H,I,L,M,O,P,Q,R,S,T,U -i "$infile_new""unique" -o itsx_sh_out_o --cpu 8 --save_regions all --partial 50 --detailed_results T -concat T -preserve T -N 1
mv itsx_sh_out_o* ITSx_o/
popd

# parse ITSx output
python3 "$script_dir/print_out_fasta.py" "$run_id" "$region"
python3 "$script_dir/print_out_fasta_o.py" "$run_id" "$region"
cat "$user_dir/seqs_out_1.fasta" "$user_dir/seqs_out_2.fasta" > "$user_dir/seqs_out.fasta"

# Chimera filtering using uchime and vsearch tools
pushd "$user_dir"
# uchime
"$program_dir/usearch" -uchime2_ref seqs_out.fasta -db "$data_dir/sanger_refs_sh.fasta" -uchimeout uchime_out.txt -strand plus -mode high_confidence
# vsearch usearch_global
vsearch -usearch_global seqs_out.fasta -db "$data_dir/sanger_refs_sh_full.unique.fasta" -strand plus -id .75 -threads 2 -uc usearch_global.full.75.map.uc --blast6out usearch_global.full.75.blast6out.txt --output_no_hits
popd

# handle all potentially chimeric sequences from uchime and usearch_global
python3 "$script_dir/exclude_chims.py" "$run_id" "$region"

# Additional quality controls - Remove low quality sequences (too short or with too many non-IUPAC symbols)
python3 "$script_dir/exclude_non_iupac.py" "$run_id" 6

# Find best matches to user’s sequences in existing SH sequence dataset using usearch_global algorithm.
pushd "$user_dir"
vsearch -usearch_global iupac_out.fasta -db "$data_dir/sanger_refs_sh.fasta" -strand plus -id .75 -threads 2 -uc closedref.75.map.uc --alnout closedref.75.aln --blast6out closedref.75.blast6out.txt --output_no_hits
popd

python3 "$script_dir/parse_vsearch_results.py" "$run_id"

# HITS: create compound clusters - Create compound clusters with user's sequences added to the existing data.

echo "Creating compound clusters"
mkdir "$user_dir/compounds"
mkdir "$user_dir/blastclust"
mkdir "$user_dir/matches"

python3 "$script_dir/create_compound_clusters.py" "$run_id"

echo "Copying compounds to blastclust folder"
pushd "$user_dir"
for nbr1 in 0 1 2 3 4 5 6 7 8 9
    do
    for nbr2 in 0 1 2 3 4 5 6 7 8 9
        do
            if ls compounds/UCL8_0$nbr1$nbr2* 1> /dev/null 2>&1
                then
                    cp compounds/UCL8_0$nbr1$nbr2* blastclust/
            fi
        done
    done
popd

bc_dir="$user_dir/blastclust"
pushd "$bc_dir"
for filename in "$bc_dir"/*.fas
    do
        mothur "#unique.seqs(fasta=$filename)"
    done
popd

# HITS: blastclust - Run blastclust for compound clusters
echo "Running blastclust for HITS sequences (compounds)..."
perl "$script_dir/blastclust_formatter2.pl" "$run_id"

# parse blastclust output
echo "Analysing BC output ..."
echo "97"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 97
echo "975"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 975
echo "98"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 98
echo "985"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 985
echo "99"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 99
echo "995"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 995
echo "100"
perl "$script_dir/analyse_BC_output.pl" "$run_id" 100

# NOHITS: Run usearch on 4 different thresholds for those sequences that didn’t match any existing SH sequence

if [ -f "$user_dir/nohits.fasta" ]
    then
        nohits_count=$(grep -c '>' "$user_dir/nohits.fasta")

        if [ "$nohits_count" -gt 0 ]
            then
                # USEARCH run
                rm -fr "$user_dir/clusters"
                mkdir "$user_dir/clusters"
                mkdir "$user_dir/clusters/clusters"
                mkdir "$user_dir/clusters/singletons"

                # 97%
                "$program_dir/usearch" -cluster_fast "$user_dir/nohits.fasta" -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -centroids "$user_dir/centroids.fasta" -uc "$user_dir/clusters_97.uc"
                python3 "$script_dir/clusterparser_preclust1.py" "$run_id"

                # 95%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_95.fasta" -id 0.95 -gapopen 0.0/0.0E -gapext 1.0/0.5E -centroids "$user_dir/centroids.fasta" -uc "$user_dir/clusters_95.uc"
                python3 "$script_dir/clusterparser_preclust2.py" "$run_id"

                # 90%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_90.fasta" -id 0.90 -gapopen 0.0/0.0E -gapext 1.0/0.5E -centroids "$user_dir/centroids.fasta" -uc "$user_dir/clusters_90.uc"
                python3 "$script_dir/clusterparser_preclust3.py" "$run_id"

                # 80%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_80.fasta" -id 0.80 -gapopen 0.0/0.0E -gapext 1.0/0.5E -centroids "$user_dir/centroids.fasta" -uc "$user_dir/clusters_80.uc"
                python3 "$script_dir/clusterparser_preclust_final.py" "$run_id"

                # NOHITS: Run blastclust for NOHITS compound clusters
                rm -fr "$user_dir"/blastclust_1
                rm -fr "$user_dir"/blastclust_out_1

                mkdir "$user_dir"/blastclust_1

                for nbr1 in 0 1 2 3 4 5 6 7 8 9
                    do
                        if ls "$user_dir"/clusters/clusters/Cluster$nbr1 1> /dev/null 2>&1
                            then
                                mv "$user_dir"/clusters/clusters/Cluster$nbr1 "$user_dir"/blastclust_1/
                        fi
                        for nbr2 in 0 1 2 3 4 5 6 7 8 9
                            do
                                if ls "$user_dir"/clusters/clusters/Cluster$nbr1$nbr2 1> /dev/null 2>&1
                                    then
                                        mv "$user_dir"/clusters/clusters/Cluster$nbr1$nbr2* "$user_dir"/blastclust_1/
                                fi
                            done
                    done

                # do blastclust clustering
                perl "$script_dir/blastclust_formatter2_1.pl" "$run_id"

                # parse blastclust output
                echo "97 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 97
                echo "975 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 975
                echo "98 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 98
                echo "985 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 985
                echo "99 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 99
                echo "995 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 995
                echo "100 (nohits)"
                python3 $script_dir/analyse_BC_output_1.py "$run_id" 100

            else
                echo "No NOHITS sequences found."
        fi
fi

# Generate output files (based on matches_*.txt files).
# parse matches files to output information about input sequences and their belonging to SHs on different thresholds
echo "Parsing SH matches ..."
perl $script_dir/parse_matches.pl "$run_id"
python3 $script_dir/parse_matches_html.py "$run_id" 97

# parse matches_1 files (nohits) to output information about input sequences and their belonging to new SHs on different thresholds
perl $script_dir/parse_matches_1.pl "$run_id"
# ToDo
# python3 $script_dir/parse_matches_1_html.py "$run_id" 97

# create Krona chart
python3 $script_dir/shmatches2kronatext.py "$run_id" 97

# export PATH=$PATH:$PROTAX/thirdparty/krona/bin
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_97.html "$user_dir"/krona_97.txt

# zip to outdata dir
pushd "$user_dir"
zip source_"$run_id".zip matches/matches_out_*.csv matches/matches_out_*.html matches/matches_1_out_*.csv err_"$run_id".log excluded_"$run_id".txt source_"$run_id"_fastanames source_"$run_id"_names krona_97.html /sh_matching/readme.txt

mv source_"$run_id".zip "$outdata_dir"/
popd

# # clean user working dir
# if [ -d "$user_dir" ]
#   then
#       rm -fr "$user_dir"
# fi

echo "End"
