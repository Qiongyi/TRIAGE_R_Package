#!/usr/bin/perl -w
use strict;
use warnings;
use LWP::Simple;
use XML::Simple;

# cpanm XML::Simple

if(@ARGV != 2) {
    print STDERR "Usage: NCBI_summary.pl Priority_epimap_rts_goslim.xls Priority_epimap_rts_annotation.xls\n";
    exit(0);
}
my ($inf, $outf)=@ARGV;

my $count=0;
my $count_miss=0;
my $count_found=0;
my $count_yes=0;
my $count_no=0;

open(OUT, ">$outf") or die "Cannot open $outf\n";
open(IN, $inf) or die "Cannot open $inf\n";
while(<IN>){
	# remove the newline character (Carriage Return + Line Feed，CRLF, "\r\n") in windows
	s/\r\n//g;
	chomp;

	if($_=~/^Ranking/){
		print OUT "$_\tGene_ID\tNCBI_summary\n";
	}elsif($_=~/\w/){
		$count++;
		my @info=split(/\t/, $_);
		my $official_symbol=$info[4];
		my $esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=$official_symbol\[sym]+AND+Homo+sapiens\[orgn]";
		my $esearch_content = get($esearch_url);
		if(defined $esearch_content){
			$count_found++;
			my $esearch_result = XML::Simple::XMLin($esearch_content);
			my $gene_id="NA";
			if (ref($esearch_result->{IdList}->{Id}) eq 'ARRAY') {
   	 			$gene_id = $esearch_result->{IdList}->{Id}->[0];
			} else {
    			$gene_id = $esearch_result->{IdList}->{Id};
			}

			# Step 2: Use ESummary with the gene ID to get the summary
			my $esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=$gene_id";
			my $esummary_content = get($esummary_url);
			if(defined $esummary_content){
				# Parse the XML to extract summary information
				my $esummary_result = XML::Simple::XMLin($esummary_content);
				my $des = $esummary_result->{DocumentSummarySet}->{DocumentSummary}->{Description};
				if($des =~/^HASH/){
					$des = "NA";
				}
				my $summary = $esummary_result->{DocumentSummarySet}->{DocumentSummary}->{Summary};
				if ($summary =~/^HASH/) {
					$summary=$des;
					$count_no++;
				}else{
					$count_yes++;
				}
				print OUT "$_\t$gene_id\t$summary\n";
			}else{
				print OUT "$_\t$gene_id\tNA\n";
				$count_no++;
			}
		}else{
			$count_miss++;
			print OUT "$_\tNA\tNA\n";
		}
	}
}
close IN;
close OUT;

print STDERR "Total number\t$count\nNumber of genes with url content\t$count_found\nNumber of genes without url content\t$count_miss\n";
print STDERR "Number of genes with summary\t$count_yes\nNumber of genes with url but no summary\t$count_no\n";


