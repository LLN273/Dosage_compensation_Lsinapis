#!/usr/bin/perl -w

use warnings;
use strict;

# Lars Höök 2017
# Filter scaffolds based on length mapped to chromosome

# Usage: perl filter_scaffolds_by_mapping_length.pl infile
# Infile format: scaffold_nr chr_nr [number of bp mapping to each chromosome separated by comma]
# (like so: scaffold_1009 1 [2680, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] +++)


####################### Apply length filter here #################################

my $filter = 200	;

# Scaffolds with "best chromosome hit" above this value will pass filter
# Examples (based on filter = 200):
		#this one will pass: scaffold_663 1 [815, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 178, 0] ++++
		#this one will be omitted: scaffold_1047 11 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 187, 0, 0, 0, 0, 0, 184, 0, 0, 0, 0] +

##################################################################################

my $infile = $ARGV[0];

my $outfile = "mapped_scaffolds_filtered_above_".$filter.".txt";
my $logfile = "mapped_scaffolds_filtered_above_".$filter.".log";
my $totalscaffolds;
my $passedscaffolds;
my $omittedscaffolds;
my $unmappedscaffolds;

open(IN, $infile) or die "Can't open $infile\n";

while (my $line = <IN>) {

	my $string = $line ;

	$string =~ s/\Q[/ /;				#remove delimiters
	$string =~ s/\Q]/ /;
	$string =~ s/\Q+/ /g;
	$string =~ s/\Q,/ /g;

	my @col = split(/\s+/, $string);		#split line by white space

	if ($col[0] =~ /scaffold_/) {			#only process lines with scaffolds

		$totalscaffolds++ ;			#log total

	my @scaffolds = splice @col, 0, 2;

	my @hits = grep($_ > 0, @col);			#grep all hits, sort in descending order and count number of hits
	my @sorted_hits = sort { $b <=> $a } @hits ;
	my $counts = @sorted_hits ;
						
		if ($counts >= 1) {			#if there are hits, apply length filter
			
		my $best_hit = $sorted_hits[0];

			if ($best_hit > $filter) {
			
			$passedscaffolds++;		#log passed

			open(OUT, ">>$outfile");
			print OUT "$line";
		
			}

			else {$omittedscaffolds++}	#log omitted
		
		}
		
		else {$unmappedscaffolds++}		#log unmapped

	}

	else {next}
	
}

close (IN);
close (OUT);

my $date = localtime();
open (LOG, ">>$logfile");
print LOG "$date\n",
	"\n", "Filter applied: $filter\n", "\n", 
	"Total scaffolds: $totalscaffolds\n",
	"Passed scaffolds: $passedscaffolds\n",
	"Ommited scaffolds: $omittedscaffolds\n",
	"Unmapped scaffolds: $unmappedscaffolds";

close (LOG);
