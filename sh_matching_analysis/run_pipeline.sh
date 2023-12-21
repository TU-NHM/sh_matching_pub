#!/bin/bash

## running the script - ./run_pipeline.sh <run_id>

## TODO LIST:
## 1. Check out log output (some steps are not logged at all)
## 2. Check if all excluded sequences are reported in del file
## 3. improve commenting inside the code
## 4. Remove unnecessary files from output during the analysis run

## check if run id has been provided
if [ -z "$1" ]
    then
        echo "No run_id argument supplied!"
        exit
fi

## Retrieve the following parameters from positional arguments:
## - 1. Run ID
## - 2. ITS region (default, "itsfull"; alternatively, "its2")
## - 3. Flag indicating whether to include the ITSx step in the analysis (default, "yes")
## - 4. Flag indicating whether to delete the user directory upon pipeline completion (default, "yes")
## - 5. Flag indicating whether to include the vsearch 100% clustering step with 96% length coverage (default, "yes")
## - 6. Flag indicating whether to conduct the usearch 99.5% (complete-linkage) distance RepS clustering step (default, "yes")
run_id=$1
region=$2
itsx_step=$3
remove_userdir=${4:-"yes"}
include_vsearch_step=${5:-"yes"}
include_usearch_05_step=${6:-"yes"}

if [ "$region" != "its2" ] && [ "$region" != "itsfull" ]; then
  echo "Setting region to itsfull"
  region="itsfull"
fi
echo "Region - $region"

if [ "$itsx_step" != "yes" ] && [ "$itsx_step" != "no" ]; then
  echo "Setting itsx_step to yes"
  itsx_step="yes"
fi

echo "ITSx - $itsx_step"

## get working directory
pwd=$(pwd)

script_dir="/sh_matching/scripts"
data_dir="/sh_matching/data"
udb_data_dir="$pwd/data_udb"
outdata_dir="$pwd/outdata"
user_dir="$pwd/userdir/$run_id"
clusters_pre_dir="$user_dir/clusters_pre"
clusters_dir="$user_dir/clusters"
program_dir="/sh_matching/programs"
infile="source_$run_id"
infile_new=$infile"_fasta"
infile_new_w_dir="$user_dir/$infile_new"

echo "Start"

echo "Starting new analysis in $user_dir"

## rm (if exists) user working dir and create as new
if [ -d "$user_dir" ]
  then
      rm -fr "$user_dir"
fi
mkdir "$user_dir"

## error log file
err_log="$user_dir/err_$run_id.log"
touch "$err_log"

## output sequences to be excluded in various analysis steps
ex_log="$user_dir/excluded_$run_id.txt"
touch "$ex_log"

## copy indata to user's workdir
cp "$pwd/indata/$infile" "$user_dir"

## replace sequence identifiers with unique codes for the analysis
python3 "$script_dir"/replace_seq_names_w_codes.py "$run_id"

## remove duplicate sequences from user’s dataset
pushd "$user_dir"
"$program_dir/vsearch/bin/vsearch" --fastx_uniques $infile_new_w_dir --fastaout "$infile_new""unique" --uc "$infile_new""uc"
popd
python3 "$script_dir"/reformat_fastx_uniques_output.py "$run_id"

if [ "$itsx_step" == "yes" ]; then
    ## Extract ITS regions, first for fungi, and then all other groups.
    pushd "$user_dir"
    mkdir ITSx
    perl "$program_dir/ITSx/ITSx" -t F -i "$infile_new""unique" -o itsx_sh_out --cpu 8 --save_regions ITS1,5.8S,ITS2 --partial 50 --detailed_results T -concat T -preserve T -N 1 --search_eval 0.1 -E 0.1 --graphical F --complement T
    mv itsx_sh_out* ITSx/

    mkdir ITSx_o
    perl "$program_dir/ITSx/ITSx" -t A,B,C,D,E,G,H,I,L,M,O,P,Q,R,S,T,U -i "$infile_new""unique" -o itsx_sh_out_o --cpu 8 --save_regions ITS1,5.8S,ITS2 --partial 50 --detailed_results T -concat T -preserve T -N 1 --graphical F --complement T
    mv itsx_sh_out_o* ITSx_o/
    popd

    ## parse ITSx output
    python3 "$script_dir/print_out_fasta.py" "$run_id" "$region"
    python3 "$script_dir/print_out_fasta_o.py" "$run_id" "$region"
    cat "$user_dir/seqs_out_1.fasta" "$user_dir/seqs_out_2.fasta" > "$user_dir/seqs_out.fasta"
else
    pushd "$user_dir"
    cp "$infile_new""unique" "seqs_out.fasta"
    popd
fi

## Chimera filtering using vsearch
pushd "$user_dir"
## vsearch usearch_global
"$program_dir/vsearch/bin/vsearch" --usearch_global "$user_dir/seqs_out.fasta" --db "$udb_data_dir/sanger_refs_sh.udb" --strand plus --id .75 --threads 8 --uc "$user_dir/usearch_global.full.75.map.uc" --blast6out "$user_dir/usearch_global.full.75.blast6out.txt" --output_no_hits
popd

## handle all potentially chimeric sequences from usearch_global
python3 "$script_dir/exclude_chims.py" "$run_id" "$region"

## Additional quality controls - Remove low quality sequences (too short or with too many non-IUPAC symbols)
python3 "$script_dir/exclude_non_iupac.py" "$run_id" 6

## Allow query sequences vary 4% in length at 100% similarity
pushd $user_dir
if [ "$include_vsearch_step" == "yes" ]; then
    echo "Running vsearch 100% clustering"
    "$program_dir/vsearch/bin/vsearch" --cluster_fast "$user_dir/iupac_out_vsearch_96.fasta" --id 1 --iddef 2 --threads 8 --uc "$user_dir/clusters_100.uc" --centroids "$user_dir/centroids_100.fasta" --query_cov 0.96 --target_cov 0.96
else
    echo "Skipping the vsearch 100% clustering step with 96% length coverage"
    "$program_dir/vsearch/bin/vsearch" --fastx_uniques "$user_dir/iupac_out_vsearch_96.fasta" --fastaout "$user_dir/centroids_100.fasta" --uc "$user_dir/clusters_100.uc"
fi
popd

## step in here with the vsearch representatives (the sequence count diff. is 9.5% for vsearch 4%)
echo "Printing out vsearch representatives"
python3 "$script_dir/select_vsearch_reps.py" "$run_id"

if [ "$include_usearch_05_step" == "yes" ]; then
    ## NEW: preclustering steps to keep only 0.5% representatives

    ## 1. usearch 97-95-90-80% clustering
    echo "usearch clustering"
    rm -fr $clusters_pre_dir
    mkdir $clusters_pre_dir
    mkdir "$clusters_pre_dir/clusters"
    mkdir "$clusters_pre_dir/singletons"

    ## 97% pre-clustering
    "$program_dir/usearch" -cluster_fast "$user_dir/iupac_out_vsearch.fasta" -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$user_dir/clusters_97_pre.uc"
    python3 "$script_dir/clusterparser_preclust1_pre.py" "$run_id"

    ## 95% pre-clustering
    "$program_dir/usearch" -cluster_fast "$user_dir/in_95_pre.fasta" -id 0.95 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$user_dir/clusters_95_pre.uc"
    python3 "$script_dir/clusterparser_preclust2_pre.py" "$run_id"

    ## 90% pre-clustering
    "$program_dir/usearch" -cluster_fast "$user_dir/in_90_pre.fasta" -id 0.90 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$user_dir/clusters_90_pre.uc"
    python3 "$script_dir/clusterparser_preclust3_pre.py" "$run_id"

    ## 80% clustering
    "$program_dir/usearch" -cluster_fast "$user_dir/in_80_pre.fasta" -id 0.80 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$user_dir/clusters_80_pre.uc"
    python3 "$script_dir/clusterparser_preclust_final_pre.py" "$run_id"

    ## remove unneeded uc, txt, and fasta files
    rm "$user_dir/clusters_97_pre.uc"
    rm "$user_dir/clusters_95_pre.uc"
    rm "$user_dir/clusters_90_pre.uc"
    rm "$user_dir/clusters_80_pre.uc"
    rm "$user_dir/clusters_out_97_pre.txt"
    rm "$user_dir/clusters_out_95_pre.txt"
    rm "$user_dir/clusters_out_90_pre.txt"
    rm "$user_dir/clusters_out_80_pre.txt"
    rm "$user_dir/in_95_pre.fasta"
    rm "$user_dir/in_90_pre.fasta"
    rm "$user_dir/in_80_pre.fasta"

    ## create a list of useach clusters
    cd "$clusters_pre_dir/clusters/"
    find . -maxdepth 1 -type f -name "Cluster*" \
        | sed 's|^./||' \
        | sort --version-sort \
        > "$clusters_pre_dir/tmp.txt"

    ## create a list of useach singletons
    cd "$clusters_pre_dir/singletons/"
    find . -maxdepth 1 -type f -name "Singleton*" \
        | sed 's|^./||' \
        | sort --version-sort \
        > "$clusters_pre_dir/singletons.txt"
    cd $pwd

    ## write vsearch clustering duplicates into duplic_seqs.txt file
    python3 "$script_dir/usearch_parser.py" "$run_id"

    ## go through 80% uclust clusters and run 97% usearch clustering if needed (if >16000 in cluster size)
    ## 2. calculate 0.5% clusters (USEARCH calc_distmx & cluster_aggd)
    touch "$user_dir/seq_mappings.txt"
    rm -fr "$clusters_pre_dir/clusters/calc_distm_out"
    mkdir "$clusters_pre_dir/clusters/calc_distm_out"

    input_95="$clusters_pre_dir/tmp.txt"
    while IFS= read -r line
    do
        fname="$clusters_pre_dir/clusters/$line"
        if [ -f "$fname" ]
            then
                result=$(grep -c ">" "$fname")
                if (( "$result" > "16000" ))
                    then
                        echo "to be split:"$line":"$result
                        "$program_dir/usearch" -cluster_fast $fname -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$clusters_pre_dir/clusters_2_90.uc"
                        ## cluster into clusters_pre/ClusterX_folder/
                        rm -fr "$clusters_pre_dir/clusters/"$line"_folder"
                        mkdir "$clusters_pre_dir/clusters/"$line"_folder"
                        mkdir "$clusters_pre_dir/clusters/"$line"_folder/clusters"
                        mkdir "$clusters_pre_dir/clusters/"$line"_folder/singletons"
                        ## parse usearch clusters
                        python3 "$script_dir/clusterparser_usearch_90_pre.py" "$run_id" "$line"
                        mkdir "$clusters_pre_dir/clusters/"$line"_folder/calc_distm_out"
                        ## calculate usearch distance matrix
                        python3 "$script_dir/calc_distm_formatter_90_pre.py" "$run_id" "$line"
                else
                    echo "sm:"$line":"$result
                    ## calculate usearch distance matrix
                    python3 "$script_dir/calc_distm_formatter_80_pre.py" "$run_id" "$line"
                fi
        fi
    done < "$input_95"

    ## 3. take 0.5% representatives as RepS, add USEARCH singletons
    python3 "$script_dir/select_core_reps_usearch.py" "$run_id"

    ## END NEW: preclustering steps to keep only 0.5% representatives
else
    cp "$user_dir/iupac_out_vsearch.fasta" "$user_dir/core_reps_pre.fasta"
fi

## Find best matches to user’s sequences in the existing SH sequence dataset using usearch_global algorithm.
pushd "$user_dir"
"$program_dir/vsearch/bin/vsearch" --usearch_global "$user_dir/core_reps_pre.fasta" --db "$udb_data_dir/sanger_refs_sh_full.udb" --strand plus --id 0.8 --threads 8 --iddef 0 --gapopen 0I/0E --gapext 2I/1E --uc "$user_dir/closedref.80.map.uc" --maxaccepts 3 --maxrejects 0
popd

python3 "$script_dir/parse_usearch_results.py" "$run_id"

## HITS: Create compound clusters (!of core dataset only!) with user's sequences added to the existing data.

echo "Creating compound clusters"
mkdir "$user_dir/compounds"
mkdir "$user_dir/matches"

python3 "$script_dir/create_compound_clusters.py" "$run_id"

## go through compound clusters and run 97% usearch clustering if needed (if >16000 in cluster size) -> calc 3.0% distance matrix to form SHs based on these
rm -fr "$user_dir/compounds/calc_distm_out"
mkdir "$user_dir/compounds/calc_distm_out"

pushd "$user_dir/compounds/"
touch "$user_dir/compounds/tmp.txt"
for nbr1 in {0..9}
    do
        for nbr2 in {0..9}
            do
                if ls UCL9_0$nbr1$nbr2* 1> /dev/null 2>&1
                    then
                        ls UCL9_0$nbr1$nbr2* >> "$user_dir/compounds/tmp.txt"
                fi
            done
    done
popd

input_95="$user_dir/compounds/tmp.txt"
while IFS= read -r line
do
    fname=$user_dir"/compounds/"$line
    if [ -f "$fname" ]
    then
        result=$(grep -c ">" "$fname")
        if (( "$result" > "16000" ))
            then
            echo "to be split:"$line":"$result
            "$program_dir/usearch" -cluster_fast $fname -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$user_dir/compounds/clusters_2_90.uc"
            ## cluster into usearch_sh/ClusterX_folder/
            rm -fr $user_dir"/compounds/"$line"_folder"
            mkdir $user_dir"/compounds/"$line"_folder"
            mkdir $user_dir"/compounds/"$line"_folder/clusters"
            mkdir $user_dir"/compounds/"$line"_folder/singletons"
            ## parse usearch clusters
            python3 $script_dir"/clusterparser_usearch_90.py" "$run_id" "$line"
            ## Calculate SHs (max 3.0% distance)
            mkdir $user_dir"/compounds/"$line"_folder/calc_distm_out"
            # calculate usearch distance matrix and generate (SH) clusters
            python3 $script_dir"/calc_distm_formatter_90.py" "$run_id" "$line"
        else
            echo "sm:"$line":"$result
            # calculate usearch distance matrix and generate (SH) clusters
            python3 $script_dir"/calc_distm_formatter_80.py" "$run_id" "$line"
        fi
  fi
done < "$input_95"

## parse usearch output
echo "Parsing usearch output ..."
echo "03"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 03
echo "025"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 025
echo "02"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 02
echo "015"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 015
echo "01"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 01
echo "005"
python3 "$script_dir/analyse_usearch_output.py" "$run_id" 005

## parse matches files to output information about input sequences and their belonging to SHs on different thresholds
echo "Parsing SH matches ..."
perl $script_dir/parse_matches.pl "$run_id"

## NOHITS: Run usearch on 4 different thresholds for those sequences that didn’t match any existing SH sequence

if [ -f "$user_dir/nohits.fasta" ]
    then
        nohits_count=$(grep -c '>' "$user_dir/nohits.fasta")

        if [ "$nohits_count" -gt 0 ]
            then
                ## USEARCH run
                rm -fr "$clusters_dir"
                mkdir "$clusters_dir"
                mkdir "$clusters_dir/clusters"
                mkdir "$clusters_dir/singletons"

                ## 97%
                "$program_dir/usearch" -cluster_fast "$user_dir/nohits.fasta" -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -uc "$user_dir/clusters_97.uc"
                python3 "$script_dir/clusterparser_preclust1.py" "$run_id"

                ## 95%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_95.fasta" -id 0.95 -gapopen 0.0/0.0E -gapext 1.0/0.5E -uc "$user_dir/clusters_95.uc"
                python3 "$script_dir/clusterparser_preclust2.py" "$run_id"

                ## 90%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_90.fasta" -id 0.90 -gapopen 0.0/0.0E -gapext 1.0/0.5E -uc "$user_dir/clusters_90.uc"
                python3 "$script_dir/clusterparser_preclust3.py" "$run_id"

                ## 80%
                "$program_dir/usearch" -cluster_fast "$user_dir/in_80.fasta" -id 0.80 -gapopen 0.0/0.0E -gapext 1.0/0.5E -uc "$user_dir/clusters_80.uc"
                python3 "$script_dir/clusterparser_preclust_final.py" "$run_id"

                ## NOHITS: Run usearch "calc_distmx & cluster_aggd" for NOHITS compound clusters
                ## go through compound clusters and run 97% usearch clustering if needed (if >16000 in cluster size) -> calc 3.0% distance matrix to form SHs based on these
                rm -fr "$clusters_dir/clusters/calc_distm_out"
                mkdir "$clusters_dir/clusters/calc_distm_out"

                pushd "$clusters_dir/clusters/"
                touch "$clusters_dir/tmp.txt"
                for nbr1 in {0..9}
                    do
                        if ls Cluster$nbr1* 1> /dev/null 2>&1
                            then
                                ls Cluster$nbr1* >> "$clusters_dir/tmp.txt"
                        fi
                    done
                    for nbr2 in {0..9}
                        do
                            if ls Cluster$nbr1$nbr2* 1> /dev/null 2>&1
                                then
                                    ls Cluster$nbr1$nbr2* >> "$clusters_dir/tmp.txt"
                            fi
                        done
                popd

                input_95="$clusters_dir/tmp.txt"
                while IFS= read -r line
                do
                    fname="$clusters_dir/clusters/$line"
                    if [ -f "$fname" ]
                    then
                        result=$(grep -c ">" "$fname")
                        if (( "$result" > "16000" ))
                            then
                            echo "to be split:"$line":"$result
                            "$program_dir/usearch" -cluster_fast $fname -id 0.97 -gapopen 0.0/0.0E -gapext 1.0/0.5E -sort other -uc "$clusters_dir/clusters_2_90.uc"
                            ## cluster into clusters/ClusterX_folder/
                            rm -fr $clusters_dir"/clusters/"$line"_folder"
                            mkdir $clusters_dir"/clusters/"$line"_folder"
                            mkdir $clusters_dir"/clusters/"$line"_folder/clusters"
                            mkdir $clusters_dir"/clusters/"$line"_folder/singletons"
                            ## parse usearch clusters
                            python3 $script_dir"/clusterparser_usearch_90_nohit.py" "$run_id" "$line"
                            ## Calculate SHs (max 3.0% distance)
                            mkdir $clusters_dir"/clusters/"$line"_folder/calc_distm_out"
                            # calculate usearch distance matrix and generate (SH) clusters
                            python3 $script_dir"/calc_distm_formatter_90_nohit.py" "$run_id" "$line"
                        else
                            echo "sm:"$line":"$result
                            # calculate usearch distance matrix and generate (SH) clusters
                            python3 $script_dir"/calc_distm_formatter_80_nohit.py" "$run_id" "$line"
                        fi
                  fi
                done < "$input_95"

                ## parse blastclust output
                echo "03 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 03
                echo "025 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 025
                echo "02 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 02
                echo "015 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 015
                echo "01 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 01
                echo "005 (nohits)"
                python3 $script_dir/analyse_usearch_output_1.py "$run_id" 005

                ## parse matches_1 files (nohits) to output information about input sequences and their belonging to new SHs on different thresholds
                perl $script_dir/parse_matches_1.pl "$run_id"
            else
                echo "No NOHITS sequences found."
        fi
fi

## merge parse_matches*.pl output into one CSV file
python3 $script_dir/merge_matches.py "$run_id"

## parse matches for html output
python3 $script_dir/parse_matches_html.py "$run_id" 005
python3 $script_dir/parse_matches_html.py "$run_id" 01
python3 $script_dir/parse_matches_html.py "$run_id" 015
python3 $script_dir/parse_matches_html.py "$run_id" 02
python3 $script_dir/parse_matches_html.py "$run_id" 025
python3 $script_dir/parse_matches_html.py "$run_id" 03

## create Krona chart
python3 $script_dir/shmatches2kronatext.py "$run_id" 005
python3 $script_dir/shmatches2kronatext.py "$run_id" 01
python3 $script_dir/shmatches2kronatext.py "$run_id" 015
python3 $script_dir/shmatches2kronatext.py "$run_id" 02
python3 $script_dir/shmatches2kronatext.py "$run_id" 025
python3 $script_dir/shmatches2kronatext.py "$run_id" 03

## export PATH=$PATH:$PROTAX/thirdparty/krona/bin
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_005.html "$user_dir"/krona_005.txt
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_01.html "$user_dir"/krona_01.txt
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_015.html "$user_dir"/krona_015.txt
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_02.html "$user_dir"/krona_02.txt
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_025.html "$user_dir"/krona_025.txt
$program_dir/krona/bin/ktImportText -o "$user_dir"/krona_03.html "$user_dir"/krona_03.txt

## zip to outdata dir
pushd "$user_dir"
zip source_"$run_id".zip matches/matches_out_*.csv matches/matches_out_*.html matches/matches_1_out_*.csv err_"$run_id".log excluded_"$run_id".txt source_"$run_id"_fastanames source_"$run_id"_names krona_*.html closedref.80.map.uc /sh_matching/readme.txt

mv source_"$run_id".zip "$outdata_dir"/
popd

## clean user working dir
if [ "$remove_userdir" == "yes" -a -d "$user_dir" ]; then
    echo "Removing user directory..."
    rm -fr "$user_dir"
else
    echo "Directory removal skipped"
fi

echo "End"
