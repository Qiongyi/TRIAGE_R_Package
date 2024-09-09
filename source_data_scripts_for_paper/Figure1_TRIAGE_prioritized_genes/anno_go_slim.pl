#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 5) {
    print STDERR "Usage: anno_go_slim.pl anno_go_slim.xls anno_description.xls gene_annotation_GeneCards_manual.txt Priority_epimap_rts_alias_DisGeNET_refseq.xls Priority_epimap_rts_annotation.xls\n";
    exit(0);
}
my ($go, $dd, $genecards, $inf, $outf)=@ARGV;

my %go;
my %des;
my %test;
open(IN, $go) or die "Cannot open $go\n";
while(<IN>){
	# remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
	s/\r\n//g;
	# remove ^M - part of the newline character (Carriage Return，CR, "\r") in windows
	#s/\r//g;
	chomp;
	next if $_=~/^hgnc_symbol/;
	my @info=split(/\t/, $_);
	$go{$info[0]}.="$info[1]/";
	$des{$info[0]}.="$info[2]/";
	$test{$info[0]}=0;
}
close IN;

# gene description from Biomart
my %dd;
open(IN, $dd) or die "Cannot open $dd\n";
while(<IN>){
	# remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
	s/\r\n//g;
	# remove ^M - part of the newline character (Carriage Return，CR, "\r") in windows
	#s/\r//g;
	chomp;
	next if $_=~/^hgnc_symbol/;
	my @info=split(/\t/, $_);
	$dd{$info[0]}=$info[1];
}
close IN;

# add gene description by manual check
open(IN, $genecards) or die "Cannot open $genecards\n";
while(<IN>){
	if($_=~/([\w\.\-]+)/){
		my $id=$1;
		my $line=<IN>;
		if($line=~/^GeneCards/){
			my $line2=<IN>;
			$line2=~s/\r\n//g;
			chomp($line2);
			$dd{$id}=$line2." [Source:GeneCards]";
		}else{
			my $line2=<IN>;
			$line=~s/\r\n//g;
			chomp($line);
			$dd{$id}=$line;
		}
	}
}
close IN;


open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	chomp;
	my @info=split(/\t/, $_);
	my $gene_description="NA";
	if(exists $dd{$info[4]}){
		$gene_description=$dd{$info[4]};
	}elsif(exists $dd{$info[1]}){
		$gene_description=$dd{$info[1]};
	}

	if($_=~/^Ranking/){
		print OUT "$_\tdescription\tgoslim_goa_accession\tgoslim_goa_description\n";
	}elsif(exists $go{$info[4]}){
		$go{$info[4]}=~s/\/$//;
		$des{$info[4]}=~s/\/$//;
		print OUT "$_\t$gene_description\t$go{$info[4]}\t$des{$info[4]}\n";
		$test{$info[4]}=1;
	}elsif(exists $go{$info[1]}){
		$go{$info[1]}=~s/\/$//;
		$des{$info[1]}=~s/\/$//;		
		print OUT "$_\t$gene_description\t$go{$info[1]}\t$des{$info[1]}\n";
		$test{$info[1]}=1;
	}else{
		#print STDERR "$_\n";
		print OUT "$_\t$gene_description\tNA\tNA\n";
	}
}
close IN;
close OUT;

foreach my $key (keys %test){
	if($test{$key}==0){
		print STDERR "$key exists in the go slim file, but not in the priority gene list. Please check.\n";
	}
}
