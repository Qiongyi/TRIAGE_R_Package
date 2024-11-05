#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 5) {
    print STDERR "Usage: anno_TF_go.pl Mus_musculus_TF.txt Mus_musculus_Cof.txt GOBP_regulation.xls inf outf\n";
    exit(0);
}
my ($tf, $cof, $go, $inf, $outf)=@ARGV;

my %tf;
open(IN, $tf) or die "Cannot open $tf: $!";
while(<IN>){
	s/\r//;
	s/\r\n//;
	my @info=split(/\t/, $_);
	next if($_=~/^Species/);
	$tf{$info[1]}=$info[3];
}
close IN;

my %cof;
open(IN, $cof) or die "Cannot open $cof: $!";
while(<IN>){
	s/\r//;
	s/\r\n//;	
	my @info=split(/\t/, $_);
	next if($_=~/^Species/);
	$cof{$info[1]}=$info[3];
}
close IN;

my %count;
my %go;
open(IN, $go) or die "Cannot open $go: $!";
while(<IN>){
	s/\r//;
	s/\r\n//;
	my @info=split(/\t/, $_);
	next if($_=~/^ID/);
	my @gene=split(/\//, $info[7]);
	foreach my $gene (@gene){
		$count{$gene}++;
		# only take the top 3 GO terms as annotation
		if($count{$gene}<=3){
			$go{$gene}.="$info[0]\($info[1]\); ";
		}
	}
}
close IN;

open(OUT, ">$outf") or die "Cannot open $outf: $!";
open(IN, $inf) or die "Cannot open $inf: $!";
while(<IN>){
	s/\r//;
	s/\r\n//;
	chomp;
	my @info=split(/[\,\t]/, $_);
	if($info[0] eq "gene_symbol"){
		print OUT "gene_symbol\tDS_value\tp_value\tannotation\n";
		next;
	}
	if($info[6]<0.01){
		my $tag=1;
		print OUT "$info[0]\t$info[5]\t$info[6]\t";
		if(exists $tf{$info[0]}){
			$tag=0;
			print OUT "TF\($tf{$info[0]}\); ";
		}
		if(exists $cof{$info[0]}){
			$tag=0;
			print OUT "Cof\($cof{$info[0]}\); ";
		}
		if(exists $go{$info[0]}){
			$tag=0;
			print OUT "$go{$info[0]}";
		}
		if($tag){
			print OUT "NA";
		}
		print OUT "\n";
	}
}
close IN;
close OUT;
