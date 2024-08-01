#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 3) {
    print STDERR "Usage: anno_gene_alias.pl gene_result.txt inf outf\n";
    exit(0);
}
my ($gene, $inf, $outf)=@ARGV;

my %hash;
my %alias;
my %loc;
open(IN, $gene) or die "Cannot open $gene\n";
while(<IN>){
	if($_=~/^tax_id/){
		next;
	}else{
		my @info=split(/\t/, $_);
		if($info[6]=~/\w/){
			$hash{$info[5]}=$info[6];
			my @ele=split(/, /, $info[6]);
			foreach my $ele (@ele){
				$alias{$ele}=$info[5];
			}
			my $id="LOC$info[2]";
			$loc{$id}="$info[5]\t$info[6]";
		}else{
			$hash{$info[5]}="NA";
			my $id="LOC$info[2]";
			$loc{$id}="$info[5]\t$info[6]";
		}
	}
}
close IN;

my $count1=0;
my $count2=0;
my $count3=0;

open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	chomp;
	my @info=split(/,/, $_);
	if($info[1] eq "Gene"){
		print OUT "Ranking\tGene\tRTS\tPriority\tOfficial_GeneSymbol\tAliases\n";
	}else{
		if(exists $hash{$info[1]}){
			$count1++;
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[1]\t$hash{$info[1]}\n";
		}elsif(exists $alias{$info[1]}){
			$count2++;
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$alias{$info[1]}\t$hash{$alias{$info[1]}}\n";
		}elsif($info[1]=~/LOC\d+/ && exists $loc{$info[1]}){
			$count2++;
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$loc{$info[1]}\n";
		}else{
			$count3++;
			print STDERR "$_\tnot found in $gene\n"; 
			print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[1]\tNA\n";
		}
	}
}
close IN;
close OUT;

print STDERR "$count1 genes match current gene symbols.\n";
print STDERR "$count2 genes match gene aliases.\n";
print STDERR "$count3 genes cannot match current gene symbols or gene aliases.\n";
