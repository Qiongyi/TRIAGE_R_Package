#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: anno_coding_noncoding_refseq.pl hg19_refseq.20240523.txt hg19.ncbiRefSeq.105.20190906.txt Priority_epimap_rts_alias_DisGeNET.xls Priority_epimap_rts_alias_DisGeNET_refseq.xls\n";
    exit(0);
}
my ($ref, $oldref, $inf, $outf)=@ARGV;

my %coding;
my %noncoding;

open(IN, $oldref) or die "Cannot open $oldref\n";
while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @info=split(/\t/, $_);
	if($info[1]=~/^NM_|^XM_|^YP_|^NP_/){
		$coding{$info[12]}=1;
	}elsif($info[1]=~/^NR_|^XR_/){
		$noncoding{$info[12]}=1;
	}else{
		print STDERR "Wrong transcript ID in $oldref: $info[1]\n"; exit;
	}
}
close IN;

open(IN, $ref) or die "Cannot open $ref\n";
while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @info=split(/\t/, $_);
	# NP and YP are used only for protein-coding genes on the mitochondrion; YP is used for human only.
	if($info[1]=~/^NM_|^XM_|^YP_|^NP_/){
		$coding{$info[12]}=1;
	}elsif($info[1]=~/^NR_|^XR_/){
		$noncoding{$info[12]}=1;
	}else{
		print STDERR "Wrong transcript ID in $ref: $info[1]\n"; exit;
	}
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	chomp;
	my @info=split(/\t/, $_);
	if($_=~/^Ranking/){
		print OUT "$_\tcoding/noncoding\n";
	}elsif(exists $coding{$info[4]} && exists $noncoding{$info[4]}){
		print OUT "$_\tboth\n";
	}elsif(exists $coding{$info[4]}){
		print OUT "$_\tcoding\n";
	}elsif(exists $noncoding{$info[4]}){
		print OUT "$_\tnoncoding\n";
	}elsif(exists $coding{$info[1]} && exists $noncoding{$info[1]}){
		print OUT "$_\tboth\n";
	}elsif(exists $coding{$info[1]}){
		print OUT "$_\tcoding\n";
	}elsif(exists $noncoding{$info[1]}){
		print OUT "$_\tnoncoding\n";
	}elsif($info[1] eq "C2orf70"){
		### manually checked: 
		# https://www.ncbi.nlm.nih.gov/gene/339778
		# CIMIP2C (C2orf70) ciliary microtubule inner protein 2C [ Homo sapiens (human) ]
		print OUT "$_\tcoding\n";
	}elsif($info[1] eq "LOC400706"){
		### manually checked: 
		# https://www.ncbi.nlm.nih.gov/gene/400706
		# Gene type: ncRNA  - Gene ID: 400706, discontinued on 23-Apr-2019
		print OUT "$_\tnoncoding\n";
	}else{
		print STDERR "$_\n";
		print OUT "$_\tNA\n";
	}
}
close IN;
close OUT;
