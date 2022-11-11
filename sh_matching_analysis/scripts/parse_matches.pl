#!/usr/bin/perl
use strict;
use warnings;

# input data
my $run_id = $ARGV[0];
if ($run_id !~ m/^[0-9]{1,}$/) {
    print "Need correct run id - number!\n";
    exit;
}

my $user_dir = "userdir/$run_id";
my $matches_dir = $user_dir . "/matches";

my $sh2compound_file = "/sh_matching/data/sh2compound_mapping.txt";
my $shs_file = "/sh_matching/data/shs_out.txt";
my $compound_file = "/sh_matching/data/compounds_out.txt";
my $centroid2sh_file = "/sh_matching/data/centroid2sh_mappings.txt";

my $accno_seqs_file = $user_dir . "/" . "source_" . $run_id . "_names";
# read in duplicates from vsearch --fastx_uniques and vsearch length coverage clustering (seq_id, parent_seq_id, cluster_name)
my $duplicate_seqs_file1 = $user_dir . "/" . "duplic_seqs.txt";
# read in duplicates from 0.5% RepS clustering (RepS_id, duplicate_seq_id)
my $duplicate_seqs_file2 = $user_dir . "/" . "seq_mappings.txt";

my %threshold_hash = ("03" => "3.0", "025" => "2.5", "02" => "2.0", "015" => "1.5", "01" => "1.0", "005" => "0.5");
my %threshold_coded_hash = ("03" => "1", "025" => "4", "02" => "2", "015" => "5", "01" => "3", "005" => "6");

# get compound and SH mappings
my %sh_ucl_hash = ();
open (SH_2_COMPOUND, $sh2compound_file);
while (<SH_2_COMPOUND>) {
    chomp $_;
    my @fields = split("\t", $_);
    $sh_ucl_hash{$fields[0]} = $fields[1];
}
close SH_2_COMPOUND;

# get refs2sh mappings for other thresholds
my %seq2sh_o_hash = ();
open (SEQ_2_SH_O, $centroid2sh_file);
while (<SEQ_2_SH_O>) {
    chomp $_;
    my @fields = split("\t", $_);
    $seq2sh_o_hash{$fields[2]}{$fields[0]} = $fields[1];
}
close SEQ_2_SH_O;

# read in SH info from shs_out.txt
my %sh_taxonomy_hash = ();
open (INFILE_SHS, $shs_file);
while (<INFILE_SHS>) {
    chomp $_;
    my @fields = split("\t", $_);
    $sh_taxonomy_hash{$fields[0]} = $fields[1];
}
close INFILE_SHS;

# read in compound info from compounds_out.txt
my %ucl_taxonomy_hash = ();
open (INFILE_COMPOUNDS, $compound_file);
while (<INFILE_COMPOUNDS>) {
    chomp $_;
    my @fields = split("\t", $_);
    $ucl_taxonomy_hash{$fields[1]} = $fields[2];
}
close INFILE_COMPOUNDS;

# read in duplicate seqs 1
my %seq_duplicate_hash = ();
open (INFILE_DUPL, $duplicate_seqs_file1);
while (<INFILE_DUPL>) {
    chomp $_;
    my @fields = split("\t", $_);
    if (!defined($seq_duplicate_hash{$fields[1]})) {
        $seq_duplicate_hash{$fields[1]} = $fields[0];
    } else {
        $seq_duplicate_hash{$fields[1]} = $seq_duplicate_hash{$fields[1]} . "," . $fields[0];
    }
}
close INFILE_DUPL;

# read in duplicate seqs 2
open (INFILE_DUPL, $duplicate_seqs_file2);
while (<INFILE_DUPL>) {
    chomp $_;
    my @fields = split(",", $_);
    if (!defined($seq_duplicate_hash{$fields[0]})) {
        if (defined($seq_duplicate_hash{$fields[1]})) {
            $seq_duplicate_hash{$fields[0]} = $fields[1] . "," . $seq_duplicate_hash{$fields[1]};
            delete($seq_duplicate_hash{$fields[1]});
        } else {
            $seq_duplicate_hash{$fields[0]} = $fields[1];
        }
    } else {
        if (defined($seq_duplicate_hash{$fields[1]})) {
            $seq_duplicate_hash{$fields[0]} = $seq_duplicate_hash{$fields[0]} . "," . $fields[1] . "," . $seq_duplicate_hash{$fields[1]};
            delete($seq_duplicate_hash{$fields[1]});
        } else {
            $seq_duplicate_hash{$fields[0]} = $seq_duplicate_hash{$fields[0]} . "," . $fields[1];
        }
    }
}
close INFILE_DUPL;

# read in sequence metadata from csv file
my %seq_id_hash = ();
open (INFILE_ACCNOS, $accno_seqs_file);
while (<INFILE_ACCNOS>) {
    chomp $_;
    my @fields = split("\t", $_);
    my $seq_id = $fields[1];
    $seq_id =~ s/i//g;
    $seq_id_hash{$seq_id} = $fields[0];
}
close INFILE_ACCNOS;

# parse files for each threshold
my $new_sh_counter = 0;
my $new_singleton_counter = 0;

my @thresholds = ("03", "025", "02", "015", "01", "005");
foreach my $threshold (@thresholds) {
    my $matches_file = $matches_dir . "/" . "matches_" . $threshold . ".txt";
    my $matches_outfile = $matches_dir . "/" . "matches_out_" . $threshold . ".csv";
    my $header = "seq_id_tmp\tseq_accno\tstatus (" . $threshold_hash{$threshold} . ")\tSH code (" . $threshold_hash{$threshold} . ")\tSH/compound taxonomy (" . $threshold_hash{$threshold} . ")\tcompound_cl_code (" . $threshold_hash{$threshold} . ")\tCompound taxonomy (" . $threshold_hash{$threshold} . ")\tduplicate_of_seq_id_tmp\tduplicate_of_seq_accno\n";

    my $new_sh_counter_th = 0;
    my $new_sh_seq_counter_th = 0;
    my $new_singleton_counter_th = 0;
    my $present_counter_th = 0;
    my %new_sh_hash_th = ();

    # open matches file
    open (MATCHES, $matches_file);
    open (MATCHES_OUT, ">", $matches_outfile);

    # print header
    print MATCHES_OUT $header;

    while (<MATCHES>) {
        chomp $_;
        my @fields = split("\t", $_);
        my $seq_id = $fields[0];
        $seq_id =~ s/i//g;
        my $best_match_seq_id = $fields[1];
        $best_match_seq_id =~ s/i//g;

        print MATCHES_OUT $seq_id . "\t" . $seq_id_hash{$seq_id} . "\t";

        my $duplicate_status = "";
        my $duplicate_sh_code = "";
        my $duplicate_ucl_code = "";
        my $duplicate_taxonomy = "";

        if ($fields[2] eq "present") {
            $present_counter_th++;
            # get SH name and taxon name for best match id
            my $sh_match_code = $seq2sh_o_hash{$threshold_coded_hash{$threshold}}{$best_match_seq_id};
            print MATCHES_OUT "present_in\t" . $sh_match_code . "\t";
            if (defined($sh_taxonomy_hash{$sh_match_code})) {
                print MATCHES_OUT $sh_taxonomy_hash{$sh_match_code} . "\t";
                $duplicate_taxonomy = $sh_taxonomy_hash{$sh_match_code};
            } else {
                print MATCHES_OUT "\t";  # SH taxonomy
            }
            $duplicate_sh_code = $sh_match_code;
            $duplicate_status = "present_in";
        } elsif ($fields[2] eq "new cluster") {
            $new_sh_seq_counter_th++;
            my $new_sh_code = "";
            if (!defined($new_sh_hash_th{$fields[4]})) {
                $new_sh_counter++;
                $new_sh_counter_th++;
                $new_sh_hash_th{$fields[4]} = $new_sh_counter;
                $new_sh_code = $new_sh_counter;
            } else {
                $new_sh_code = $new_sh_hash_th{$fields[4]};
            }
            print MATCHES_OUT "new_sh_in\t" . $new_sh_code . "\t";
            # SH taxonomy
            if (defined($fields[3]) && defined($ucl_taxonomy_hash{$fields[3]})) {
                print MATCHES_OUT $ucl_taxonomy_hash{$fields[3]} . "\t";
                $duplicate_taxonomy = $ucl_taxonomy_hash{$fields[3]};
            } else {
                print MATCHES_OUT "\t";
            }
            $duplicate_sh_code = $new_sh_code;
            $duplicate_status = "new_sh_in";
        } elsif ($fields[2] eq "singleton") {
            $new_singleton_counter++;
            $new_singleton_counter_th++;
            print MATCHES_OUT "new_singleton_in\ts". $new_singleton_counter . "\t";
            # SH taxonomy
            if (defined($fields[3]) && defined($ucl_taxonomy_hash{$fields[3]})) {
                print MATCHES_OUT $ucl_taxonomy_hash{$fields[3]} . "\t";
                $duplicate_taxonomy = $ucl_taxonomy_hash{$fields[3]};
            } else {
                print MATCHES_OUT "\t";
            }
            $duplicate_sh_code = "s" . $new_singleton_counter;
            $duplicate_status = "new_singleton_in";
        }
        if (!defined($fields[3])) {
            print "Not defined - $seq_id\n";
            print MATCHES_OUT "\t\t";
        } else {
            print MATCHES_OUT $fields[3] . "\t";
            $duplicate_ucl_code = $fields[3];
            # UCL taxonomy
            if (defined($ucl_taxonomy_hash{$fields[3]})) {
                print MATCHES_OUT $ucl_taxonomy_hash{$fields[3]} . "\t";
            } else {
                print MATCHES_OUT "\t";
            }
        }
        print MATCHES_OUT "\t\t\n";
        # check duplicates
        if (defined($seq_duplicate_hash{$fields[0]})) {
            my @fields2 = split(",", $seq_duplicate_hash{$fields[0]});
            for (my $k=0; $k<scalar(@fields2); $k++) {
                my $dupl = $fields2[$k];
                $dupl =~ s/i//g;
                # print each duplicate in a separate row
                print $dupl . "\t" . $seq_id_hash{$dupl} . "\tDuplicate of " . $seq_id_hash{$seq_id} . "\n";
                print MATCHES_OUT $dupl . "\t" . $seq_id_hash{$dupl} . "\t";
                print MATCHES_OUT $duplicate_status . "\t" . $duplicate_sh_code . "\t" . $duplicate_taxonomy . "\t" . $duplicate_ucl_code . "\t" . $ucl_taxonomy_hash{$duplicate_ucl_code} . "\t" . $seq_id . "\t" . $seq_id_hash{$seq_id} . "\n";
            }
        }
    }
    close MATCHES;
    close MATCHES_OUT;

    print "\nNo. of sequences present in SHs ($threshold): " . $present_counter_th . "\n";
    print "No. of sequences in new SHs ($threshold): " . $new_sh_seq_counter_th . "\n";
    print "No. of new SHs ($threshold): " . $new_sh_counter_th . "\n";
    print "No. of new singleton SHs ($threshold): " . $new_singleton_counter_th . "\n";
}

print "\nTotal no. of new SHs: " . $new_sh_counter . "\n";
print "Total no. of new singleton SHs: " . $new_singleton_counter . "\n\n";

exit;
