#!/usr/bin/perl -w

use warnings;
use strict;

# Use assigned as A or Z as infile
# Sort genes by highest FPKM value (male or female) for each gene and print to new file
# The outfile can then be used to split data into quartiles while keeping the FPKM values for each gene and sex separated

# Lars Höök 2017


if (@ARGV != 1) {

print "Usage: perl quartile_analysis_1.pl in.file\n";
exit;

} 

my ($infile) = @ARGV;
my $outfile = "sorted-$infile";

open (IN, $infile) or die "Can't open $infile\n";

my %data;



while (my $line = <IN>) {
	
	chomp $line;
	my @col = split(/\t/, $line);
	my $FPKM1 = $col[3];
	my $FPKM2 = $col[4];
	my $chromosome = $col[5];

	if ($col[0] eq "gene_id") {			#print headers
	
		open (OUT, ">>$outfile");
		print OUT "gene_id	$FPKM1	$FPKM2	max_FPKM\n";

		}

		elsif ($chromosome eq "A") {next}		
		elsif ($chromosome eq "Z") {			#select Z-linked genes

		splice @col,1,2;		
		my $Z = pop @col;

		my $key = join '	', @col;		#make hash key array
		
		my @sorted_FPKM = sort { $b <=> $a } @col[1,2];	#sort the two FPKM values based on highest of male or female
				
		$data{$key} = $sorted_FPKM[0];			#set highest FPKM as hash value for sorting

	}

}

	open (OUT, ">>$outfile");

	foreach my $key (sort { $data{$a} <=> $data{$b} } keys %data) { 	#sort all genes by fpkm
	
	print OUT "$key	$data{$key}\n";

}


