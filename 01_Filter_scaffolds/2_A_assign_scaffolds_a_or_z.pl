#!/usr/bin/perl -w

use strict;
use warnings;

# Lars Höök 2017
# Assign scaffolds by fraction of total mapped length belonging to either A or Z

# Useage: perl 2_A_assign_scaffolds_a_or_z.pl infile
# Infile format: scaffold_nr chr_nr [number of bp mapping to each chromosome separated by comma]
# Z chromosome should be last in array (else switch pop to shift at lines 59-60)
# (like so: scaffold_1009 1 [2680, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] +++)


####################### Apply fraction filter here ###############################

my $filter = 0.95	;

# Scaffolds pass if the mapped length is autosomal or Z-linked above (or equal to) the given filter  
# Examples (based on filter = 0.95):
		#this one will be assigned A: scaffold_128 1 [7122, 0, 0, 0, 0, 0, 0, 162, 0, 325, 0, 0, 163, 330, 333, 0, 0, 501, 0, 167, 0]  A
		#this one will be assigned Z: scaffold_0 21 [0, 74, 146, 0, 0, 0, 228, 0, 165, 71, 0, 270, 204, 0, 0, 71, 0, 70, 0, 0, 53844]  Z 
		#this one will be omitted: scaffold_766 21 [0, 0, 0, 0, 0, 366, 0, 0, 250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 481] ++++

##################################################################################

my ($infile) = $ARGV[0] ;

my $prefix = $infile;
	$prefix =~ s/\.txt//;

my $outfile = $prefix."_assigned_a_or_z_by_".$filter.".txt";
my $logfile = $prefix."_assigned_a_or_z_by_".$filter.".log";
my $totalscaffolds;
my $assignedscaffolds;
my $A_scaffolds;
my $Z_scaffolds;
my $omittedscaffolds;

open(IN, $infile) or die "Can't open $infile\n";

while (my $line = <IN>) {

	$totalscaffolds++ ;
	$line =~ s/\Q+//g;
	chomp($line);	
	
	my $string = $line ;
						
	$string =~ s/\Q[/ /;				#clean up string
        $string =~ s/\Q]/ /;
	$string =~ s/\Q,/ /g;
						
	my @col = split(/\s+/, $string);		#split by white space and remove first two columns, leaving only chromosome mapping
	my @scaffold = splice @col, 0, 2;
				
	my $sum;					#sum values of all chromosomes for calculation of fraction
	grep { $sum += $_ } @col;
		
	my $z_chr = pop(@col);				#assign the last value in @col as Z-chromosome, rest as autosomes
	#my $z_chr = shift(@col);			#if Z chromosome is first in array, switch from pop to shift here


		if ($z_chr == 0) {			#if Z = 0 scaffold is A
				
		$A_scaffolds++;
		$assignedscaffolds++;			#log assigned

		open(OUT, ">>$outfile");
		print OUT "$line"." A\n";
				
		}
			
		else { 					#else it can be Z or A
								
		my $fraction = ($z_chr / $sum);

			if ($fraction >= $filter) {	#if fraction of Z is above filter, scaffold is Z
	
			$Z_scaffolds++;			#log Z
			$assignedscaffolds++;		#log assigned	

			open(OUT, ">>$outfile");
			print OUT "$line"." Z\n";
				
			}						

			elsif ($fraction <= 1-$filter) {#if fraction of Z is below 1-filter, scaffold is A

			$A_scaffolds++;			#log A
			$assignedscaffolds++;		#log assigned
					
			open(OUT, ">>$outfile");
			print OUT "$line"." A\n";
														
			}
					
			else{$omittedscaffolds++}

		}		
			
}


close (IN);
close (OUT);

my $date = localtime();
open (LOG, ">>$logfile");
print LOG "$date\n",
	"\n", "Filter applied: $filter\n", "\n", 
	"Total scaffolds: $totalscaffolds\n",
	"Assigned scaffolds: $assignedscaffolds\n",
	"Autosomal scaffolds: $A_scaffolds\n",
	"Z-linked scaffolds: $Z_scaffolds\n",
	"Ommited scaffolds: $omittedscaffolds\n";

close (LOG);
