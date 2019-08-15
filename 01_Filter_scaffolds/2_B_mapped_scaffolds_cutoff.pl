#!/usr/bin/perl -w

use strict;
use warnings;

# Lars Höök 2017
# Assign scaffolds to individual chromosomes

# Useage: perl assign_scaffolds_to_chromosomes.pl infile
# Infile format: scaffold_nr chr_nr [number of bp mapping to each chromosome separated by comma]
# (like so: scaffold_1009 1 [2680, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] +++)


####################### Apply filter and Z chromosome nr here ###############################

my $filter = 0.9	;
my $z_chromosome = 21	;

# Scaffolds are assigned if a fraction (given by the filter) of the mapped length belongs to one chromosome
# Examples (based on filter = 0.9):
		#this one will be assigned to chromosome 1: scaffold_210 1 [15611, 0, 147, 85, 80, 78, 0, 0, 0, 0, 0, 80, 0, 0, 0, 0, 0, 0, 204, 0, 0]
		#this one will be omitted: scaffold_244 1 [1751, 0, 0, 146, 0, 0, 0, 0, 94, 822, 0, 72, 72, 82, 0, 108, 0, 74, 0, 0, 87]

# Z chromosome will be put first (for plotting): Z, 1, 2, 3...N

#############################################################################################


my ($infile) = $ARGV[0] ;

my $prefix = $infile;
$prefix =~ s/\.txt//;

my $outfile = $prefix."_assigned_to_chromosome_by_".$filter.".txt";
my $logfile = $prefix."_assigned_to_chromosome_by_".$filter.".log";
my $totalscaffolds;
my $assignedscaffolds;
my $omittedscaffolds;

open(IN, $infile) or die "Can't open $infile\n";

while (my $line = <IN>) {

	$totalscaffolds++ ;
	$line =~ s/\Q+//g;
	chomp($line);
	
	my @col_line = split(/\s+/, $line);
	
		if($col_line[1] == $z_chromosome) {	#set Z chromosome number to 0 for placing Z first when plotting
				
		$col_line[1] = "0";
		
		}
		
		elsif($col_line[1] < 10) {		#add 0 before <10 chromosome numbers (if needed) to ensure correct sorting
	
		$col_line[1] = "0" . $col_line[1];
	
		}	
	
	my $string = $line ;
						
	$string =~ s/\Q[/ /;				#clean up string
        $string =~ s/\Q]/ /;
	$string =~ s/\Q,/ /g;
						
	my @col_string = split(/\s+/, $string);		#split by white space and remove first two columns, leaving only chromosome mapping
	my @scaffold = splice @col_string, 0, 2;
						
	my @hits = grep( $_ > 0, @col_string );		#grep all hits and sort in descending order
	my @sorted_hits = sort { $b <=> $a } @hits ;	

	my $best_hit = $sorted_hits[0];			#store largest value as best hit
	my $sum;
	grep { $sum += $_ } @sorted_hits;		#sum all values
	my $fraction = ($best_hit / $sum);		#divide largest value by sum
			
		if ($fraction >= $filter) {		#assign scaffold if fraction >= filter
				
		$assignedscaffolds++;			#log assigned

		open(OUT, ">>$outfile");
		print OUT "@col_line\n";					
			
		}

		else {$omittedscaffolds++}		#log omitted

}

close (IN);
close (OUT);


my $date = localtime();

open (LOG, ">>$logfile");
print LOG "$date\n",
	"\n", "Filter applied: $filter\n", "\n", 
	"Total scaffolds: $totalscaffolds\n",
	"Assigned scaffolds: $assignedscaffolds\n",
	"Ommited scaffolds: $omittedscaffolds\n";

close (LOG);
