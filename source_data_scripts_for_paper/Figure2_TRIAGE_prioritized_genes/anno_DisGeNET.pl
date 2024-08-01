#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 3) {
    print STDERR "Usage: anno_DisGeNET.pl gene_associations.tsv Priority_epimap_rts_alias.xls Priority_epimap_rts_alias_DisGeNET.xls\n";
    exit(0);
}
my ($dis, $inf, $outf)=@ARGV;

my %hash;
open(IN, $dis) or die "Cannot open $dis\n";
while(<IN>){
	chomp;
	s/^\s+//;
	my @info=split(/\t/, $_);
	if($_=~/^geneId/){
		next;
	}else{
		my @info=split(/\t/, $_);
		if($info[5]=~/\w/){
			$hash{$info[1]}="$info[5]\t$info[7]";
		}else{
			$hash{$info[1]}="NA\t$info[7]";
		}
	}
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	chomp;
	my @info=split(/\t/, $_);
	if($_=~/^Ranking/){
		print OUT "$_\tprotein_class_name\tNofDisease\n";
	}elsif(exists $hash{$info[4]}){
		print OUT "$_\t$hash{$info[4]}\n";
	}else{
		print OUT "$_\tNA\tNA\n";
	}
}
close IN;
close OUT;
