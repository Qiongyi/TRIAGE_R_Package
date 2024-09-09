#!/usr/bin/perl -w
use strict;
use warnings;

if (@ARGV != 3) {
    die "Usage: anno_FUMA.pl FUMA_DIR inf outf\n";
}

my ($dir, $inf, $outf) = @ARGV;

my $tf_target="$dir/TF_targets.gmt";
my $gwas_catalog="$dir/GWAScatalog.gmt";

my %tf;
open(IN, $tf_target) or die "Cannot open $tf_target\n";
while(<IN>){
    # remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
    s/\r\n//g;
    chomp;    
    my @info=split(/\t/, $_);
    for(my $i=2; $i<=$#info; $i++){
        $tf{$info[$i]}.="$info[0]/";
    }
}
close IN;

my %gwas;
open(IN, $gwas_catalog) or die "Cannot open $gwas_catalog\n";
while(<IN>){
    # remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
    s/\r\n//g;
    chomp;
    my @info=split(/\t/, $_);
    for(my $i=2; $i<=$#info; $i++){
        $gwas{$info[$i]}.="$info[0]/";
    }
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
    # remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
    s/\r\n//g;
    chomp;
    my @info=split(/\t/, $_);
    my $tmp=join"\t", @info[0..$#info-1];
    if($_=~/^Ranking/){
        print OUT "$tmp\tTF_target\tGWAS_catalog\tNCBI_summary\n";
    }else{
        if(!exists $tf{$info[$#info-1]}){
            $tf{$info[$#info-1]}="NA";
        }   
        if(!exists $gwas{$info[$#info-1]}){
            $gwas{$info[$#info-1]}="NA";
        }
        $tf{$info[$#info-1]}=~s/\/$//g;
        $gwas{$info[$#info-1]}=~s/\/$//g;
        print OUT "$tmp\t$tf{$info[$#info-1]}\t$gwas{$info[$#info-1]}\t$info[$#info]\n";
    }
}
close IN;
close OUT;
