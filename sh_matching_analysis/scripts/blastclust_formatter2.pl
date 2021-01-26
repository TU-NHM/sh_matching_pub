#!/usr/bin/perl -w
use strict;

# input data
my $run_id = $ARGV[0];
if ($run_id !~ m/^[0-9]{1,}$/) {
    print "Need correct run id - number!\n";
    exit;
}
my $user_dir = "userdir/$run_id";
my $folder = $user_dir . "/blastclust";
my $bc_location = "/sh_matching/programs/blast-2.2.26/bin/blastclust";

my $check_folder = $folder . "/*.unique.fas";
my @files = glob($check_folder);

for (my $k=0; $k<scalar(@files); $k++) {
    my $name = $files[$k];
    my $ucl_code = (split '\/', $name)[-1];
    my $name_out_97 = $name . "_out_97";
    my $name_out_975 = $name . "_out_975";
    my $name_out_98 = $name . "_out_98";
    my $name_out_985 = $name . "_out_985";
    my $name_out_99 = $name . "_out_99";
    my $name_out_995 = $name . "_out_995";
    my $name_out_100 = $name . "_out_100";

    my $system_command1 = $bc_location . " -i " . $name . " -S 97 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_97 . " -p F";
    system($system_command1);
    my $system_command15 = $bc_location . " -i " . $name . " -S 97.5 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_975 . " -p F";
    system($system_command15);
    my $system_command2 = $bc_location . " -i " . $name . " -S 98 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_98 . " -p F";
    system($system_command2);
    my $system_command25 = $bc_location . " -i " . $name . " -S 98.5 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_985 . " -p F";
    system($system_command25);
    my $system_command3 = $bc_location . " -i " . $name . " -S 99 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_99 . " -p F";
    system($system_command3);
    my $system_command35 = $bc_location . " -i " . $name . " -S 99.5 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_995 . " -p F";
    system($system_command35);
    my $system_command4 = $bc_location . " -i " . $name . " -S 100 -L 0.95 -b F -a 8 -e F -W 16 -o " . $name_out_100 . " -p F";
    system($system_command4);
}

exit;
