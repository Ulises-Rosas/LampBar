#!/usr/bin/perl -w

use strict;
use warnings;

my $input_fasta = $ARGV[0];

unless(open(IN, "$input_fasta")){
	print "\nError upon opening \"$input_fasta\"\n\n";
	exit;
	}

my @lines  = <IN>;
close IN;

my $line =  join("", @lines);
my $results = "";

for $line (@lines){
	if($line =~ m/^>/){
	$results .=  "\n$line"; 
	}else{
	$line =~ s/\s//g;
	$results .=  $line;
	}
}

$results =~ s/^\n//g;
print "$results\n";
exit;
