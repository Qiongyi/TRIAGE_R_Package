#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 2) {
    print STDERR "Usage: Runtime_TRIAGEcluster2.pl outf outf2(failed_info)\n";
    exit(0);
}

my ($outf, $outf2) = @ARGV;

my %jobs = (
    '10G'  => [5000, 10000, 20000],
    '30G' => [30000, 40000, 50000, 60000],
    '100G' => [70000, 80000, 90000, 100000],
);

my @job_ids;
my %hash;

foreach my $memory (keys %jobs) {
    my @sample_sizes = @{$jobs{$memory}};

    foreach my $size (@sample_sizes) {
        my $script = "TRIAGEcluster_${size}_cells.R";
        my $job_output = `run_bunya_R.pl $memory 1 $script`;
        
        if ($job_output =~ /Submitted batch job (\d+)/) {
            my $job_id = $1;
            print "Submitted job ID for ${size} cells: $job_id\n";
            push @job_ids, $job_id;
            $hash{$job_id}=$size;
        }
    }
}

my $print_header = 1;
if (-e $outf) {
    my $size = -s $outf;
    if ($size > 0) {
        $print_header = 0;
    }
}

open(OUT, ">>$outf") or die "Cannot open $outf\n";
open(OUT2, ">>$outf2") or die "Cannot open $outf2\n";

if ($print_header) {
    print OUT "Cell_number\tJob_id\tMax_memory\n";
}

my %completed_jobs;
while (@job_ids) {
    foreach my $job_id (@job_ids) {
        my $state_output = `sacct -j $job_id --format=State --noheader`;
        my @states = split(/\n/, $state_output);
        #my %unique_states = map { $_ => 1 } @states;
        my %unique_states = map { s/^\s+|\s+$//g; $_ => 1 } @states;

        if (exists($unique_states{"COMPLETED"}) && !$completed_jobs{$job_id}) {
            my $seff_output = `seff $job_id`;
            if ($seff_output =~ /Memory Utilized:\s+([\d\.]+\s+\w+)/) {
                my $max_memory = $1;
                print OUT "$hash{$job_id}\t$job_id\t$max_memory\n";
            }
            $completed_jobs{$job_id} = 1;  # Mark as completed
        }
        elsif (exists($unique_states{"FAILED"}) || exists($unique_states{"CANCELLED"}) || exists($unique_states{"OUT_OF_ME+"})) {
            print OUT2 "Job $job_id ($hash{$job_id}) has failed or been cancelled.\n";
            print STDERR "Job $job_id ($hash{$job_id}) has failed or been cancelled.\n";
            $completed_jobs{$job_id} = 1;  # Mark as completed
        }
    }
    @job_ids = grep { !$completed_jobs{$_} } @job_ids;
    sleep(60); 
}

close OUT;
close OUT2;