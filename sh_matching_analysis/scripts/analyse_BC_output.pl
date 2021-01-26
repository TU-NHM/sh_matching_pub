#!/usr/bin/perl
use strict;
use warnings;

# input data
my $run_id = $ARGV[0];
my $threshold = $ARGV[1];  # 97, 975, 98, 985, 99, 995, 100
if ($run_id !~ m/^[0-9]{1,}$/) {
    print "Need correct run id - number!\n";
    exit;
}

my $user_dir = "userdir/$run_id";
my $matches_dir = $user_dir . "/matches";
my $sh2compound_file = "/sh_matching/data/sh2compound_mapping.txt";

# open log file
my $log_file = $user_dir . "/" . "err_" . $run_id . ".log";
open (LOG, ">>", $log_file);

# open matches file
my $matches_file = $matches_dir . "/" . "matches_" . $threshold . ".txt";
open (MATCHES, ">", $matches_file);

# read in results from prev version (if singleton, no need to check it again)
my $prev_file = "";
if ($threshold eq "975") {
    $prev_file = $matches_dir . "/" . "matches_97.txt";
} elsif ($threshold eq "98") {
    $prev_file = $matches_dir . "/" . "matches_975.txt";
} elsif ($threshold eq "985") {
    $prev_file = $matches_dir . "/" . "matches_98.txt";
} elsif ($threshold eq "99") {
    $prev_file = $matches_dir . "/" . "matches_985.txt";
} elsif ($threshold eq "995") {
    $prev_file = $matches_dir . "/" . "matches_99.txt";
} elsif ($threshold eq "100") {
    $prev_file = $matches_dir . "/" . "matches_995.txt";
}

my %prev_hash = ();
if ($prev_file ne "") {
    open (PREV_FILE, $prev_file);
    while (<PREV_FILE>) {
        chomp $_;
        my @fields = split("\t", $_);
        if ($fields[2] eq "singleton") {
            $prev_hash{$fields[0]} = "1";
        }
    }
    close PREV_FILE;
}

# get compound and SH mappings
my %sh_ucl_hash = ();
open (SH_2_COMPOUND, $sh2compound_file);
while (<SH_2_COMPOUND>) {
    chomp $_;
    my @fields = split("\t", $_);
    $sh_ucl_hash{$fields[0]} = $fields[1];
}
close SH_2_COMPOUND;

my @names_files = glob($user_dir . "/blastclust/*.names");

my %seq_mapping_hash = ();
my %seq_mapping_hash_op = ();
foreach my $file (@names_files) {
    open (NAMES_FILE, "$file");
    while (<NAMES_FILE>) {
        chomp $_;
        my @fields = split("\t", $_);
        if ($fields[0] ne $fields[1]) {
            $seq_mapping_hash{$fields[0]} = $fields[1];
        }
    }
    close NAMES_FILE;
}

# read seq ids and their UCL belongings into hash
my %seq_ucl_mapping_hash = ();
my @files = glob($user_dir . "/blastclust/*.unique.fas_out_" . $threshold);
for (my $k=0; $k<scalar(@files); $k++) {
    my $name = $files[$k];
    open (UCL_FILE, "$name");
    while (<UCL_FILE>) {
        chomp $_;
        my @ucl_file_fields = split(" ", $_);
        for (my $k=0; $k < scalar(@ucl_file_fields); $k++) {
            $seq_ucl_mapping_hash{$ucl_file_fields[$k]} = $name;
            if (defined($seq_mapping_hash{$ucl_file_fields[$k]})) {
                my @mappings = split(",", $seq_mapping_hash{$ucl_file_fields[$k]});
                for (my $l=1; $l < scalar(@mappings); $l++) {
                    $seq_ucl_mapping_hash{$mappings[$l]} = $name;
                    $seq_mapping_hash_op{$mappings[$l]} = $ucl_file_fields[$k];
                }
            }
        }
    }
    close UCL_FILE;
}

# read in best matches (seq and SH)
my $tmp_infile = $user_dir . "/closedref.75.map.uc";
open (INFILE, $tmp_infile);
while (<INFILE>) {
    chomp $_;
    my @fields = split("\t", $_);
    if ($fields[0] eq "H") {
        my @fields3 = split("_", $fields[9], 2);
        my $query = $fields[8];
        my $subject = $fields3[0];
        my $subject_sh = $fields3[1];
        if (!defined($prev_hash{$query})) {
            print MATCHES $query . "\t" . $subject . "\t";
            if ($fields[3] ne "100.0") {
                my $tmp_folder = $user_dir . "/blastclust/" . $sh_ucl_hash{$subject_sh} . ".unique.fas_out_" . $threshold;
                my $check = `grep '$query' $tmp_folder | grep '$subject'`;
                if (defined($seq_mapping_hash_op{$query})) {
                    $check = `grep '$seq_mapping_hash_op{$query}' $tmp_folder | grep '$subject'`;
                } elsif (defined($seq_mapping_hash_op{$subject})) {
                    $check = `grep '$query' $tmp_folder | grep '$seq_mapping_hash_op{$subject}'`;
                } elsif (defined($seq_mapping_hash_op{$query}) && defined($seq_mapping_hash_op{$subject})) {
                    $check = `grep '$seq_mapping_hash_op{$query}' $tmp_folder | grep '$seq_mapping_hash_op{$subject}'`;
                }

                if (defined($check) && ($check ne "")) {
                    print MATCHES "present\t" . $sh_ucl_hash{$subject_sh} . "\t" . "\n";
                } else {
                    my $check_cmd_1 = "grep '$query' " . $seq_ucl_mapping_hash{$query};
                    if (defined($seq_mapping_hash_op{$query})) {
                        $check_cmd_1 = "grep '$seq_mapping_hash_op{$query}' " . $seq_ucl_mapping_hash{$query};
                    }
                    $check = `$check_cmd_1`;
                    chomp $check;
                    $check =~ s/[\s\r\n]+$//;
                    if (!defined($check) || ($check eq "")) {
                        print MATCHES "missing\t\t\n";
                    } else {
                        if ($check eq $query) {
                            print MATCHES "singleton\t" . $sh_ucl_hash{$subject_sh} . "\t" . "\n";
                        } else {
                            print MATCHES "new cluster\t" . $sh_ucl_hash{$subject_sh} . "\t" . $check . "\n";
                        }
                    }
                }
            } else {
                print MATCHES "present\t" . $sh_ucl_hash{$subject_sh} . "\t" . "\n";
            }
        } else {
            print MATCHES $query . "\t" . $subject . "\t";
            print MATCHES "singleton\t" . $sh_ucl_hash{$subject_sh} . "\t" . "\n";
        }
    } else {
        print LOG "BC\tno hit for " . $fields[8] . "\n";
    }
}
close INFILE;

close LOG;

exit;
