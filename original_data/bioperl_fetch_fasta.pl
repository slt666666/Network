#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;

#--------------------------------------------------------------------------
# bioperl_fetch_fasta.pl
#--------------------------------------------------------------------------
#
#
# This program retrives fasta sequences from a fasta multiple sequence 
# file. The list of sequence names to retrieve must be in a text file with
# each name per line.
#
# Uses Bio::SeqIO from Bioperl
#
#
# Usage: perl bioperl_fetch_fasta.pl <name_list_file> <input_fasta_file_to_search> <output_file>
#    
#
# The Sainsbury Laboratory, Norwich NR4 7UH, UK
# joe.win@tsl.ac.uk
# 2010©
#
#--------------------------------------------------------------------------

my ($seqNamesFile, $fastafile, $outFile) = @ARGV;

my $usage = 'Usage: perl fetch_fasta.pl <name_list_file> <input_fasta_file_to_search> <output_file>';
unless ($seqNamesFile && $fastafile && $outFile) {die "$usage\nAll three files must be supplied!\n"}
unless (-e $seqNamesFile) {die "Can't open $seqNamesFile: $!"}
unless (-e $fastafile) {die "Can't open $fastafile: $!"}

my $in = Bio::SeqIO->new(-file => $fastafile,
                      -format => 'fasta');

my $out = Bio::SeqIO->new(-file =>">$outFile",
                        -format => 'fasta');

my %seqNames = ();
open (NAMES, $seqNamesFile) || die "Can't open $seqNamesFile: $!";

# Reads in the names of the sequences to search into an array
while (<NAMES>) {
    chomp;
    next if (/^\s*$/);
    $seqNames{$_} = $_;
}
close NAMES;

my %found_list = ();
my %not_found_list = ();

# Reads in the sequences from the fasta file to search and match names
while (my $obj = $in->next_seq) {
	my $this_id = $obj->display_id;
	if (exists $seqNames{$this_id}) {
		$out->write_seq($obj);
		$found_list{$this_id} = $this_id;
		delete ($seqNames{$this_id});
	}
}


my $not_found_count = keys %seqNames;
if ($not_found_count == 0) {
	print "All sequences were found and extracted into \"$outFile\"\n";
	exit;
}
	
print "Found the following sequences: \n";
Print_List (\%found_list);

print "Can\'t find the following sequences: \n";
Print_List (\%seqNames);

exit;


#--------------------------------------------------------------------------
# Subroutines
#--------------------------------------------------------------------------


sub Print_List {
	# prints out all values in a hash passed as a reference
	# prints 5 items per line
	my $hash_ref = shift @_;
	my $count = 0;
	foreach my $key(keys %$hash_ref) {
		print "$$hash_ref{$key}\t";
		$count++;
		if ($count == 5) {
			print "\n";
			$count = 0;
		}
	}
	print "\n";
}

#--------------------------------------------------------------------------
