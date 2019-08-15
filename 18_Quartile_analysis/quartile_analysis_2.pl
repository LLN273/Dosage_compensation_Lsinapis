#! usr/bin/perl -w

use warnings;
use strict;

# Split dataset into quartiles based on max of male or female FPKM

# Lars Höök 2017

if (@ARGV != 5) {

print "Usage: perl quartile_analysis_2.pl in.file q1_size q2_size q3_size q4_size\n";
exit;

}

my $infile = $ARGV[0];
my $outfile = "quartiles-$infile";
my $q1 = $ARGV[1];
my $q2 = $ARGV[2];
my $q3 = $ARGV[3];
my $q4 = $ARGV[4];

open (IN, $infile) or die "Can't open infile\n";

while (my $line = <IN>) {

	my @col = split(/\t/, $line);
	my $gene_id = $col[0];
	my $FPKM1 = $col[1];
	my $FPKM2 = $col[2];

	if ($gene_id eq "gene_id") {	
			
		open (OUT, ">>$outfile") ;
		print OUT "gene_id	$FPKM1	$FPKM2	quartile\n" ;		#print headers
		
		}
		
		else {							#make quartiles
		
		if ($q1 > 0) { 
			
		$q1--;
			
		open(OUT, ">>$outfile");
		print OUT "$gene_id	$FPKM1	$FPKM2	q_1\n";
			
		next;

		}

		elsif ($q2 > 0) {
			
		$q2--;

		open(OUT, ">>$outfile");
		print OUT "$gene_id	$FPKM1	$FPKM2	q_2\n";
				
		next;
		
		}
			
		elsif ($q3 > 0) {
			
		$q3--;

		open(OUT, ">>$outfile");
		print OUT "$gene_id	$FPKM1	$FPKM2	q_3\n";		
		
		next;

		}

		elsif ($q4 > 0) {
			
		$q4--;

		open(OUT, ">>$outfile");
		print OUT "$gene_id	$FPKM1	$FPKM2	q_4\n";		

		next;
			
		}
	}
}

