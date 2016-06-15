#!/usr/bin/env perl

use strict;
use warnings;
use v5.10;

use Getopt::ArgParse;

my $parser = Getopt::ArgParse->new_parser(
    prog        => 'parse-treks-output.pl',
    description => 'Parse MSA and tab-delimited data from T-Reks '
);
$parser->add_arg('input',  required => 1);
$parser->add_arg('tabout', required => 1);
$parser->add_arg('msaout', required => 1);

my $args = $parser->parse_args();

open(INPUT, "<", $parser->input)
    or die "Can't open < $parser->input: $!";

open(TABOUT, ">", "$parser->tabout")
    or die "Can't open > $parser->tabout: $!";

open(MSAOUT, ">", "$parser->msaout")
    or die "Can't open > $parser->msaout: $!";

my $seqid;
my $msa_count;
while(<INPUT>){

    # Find proteins with no repeats found
    if(/not found.*>(\S+)/) {
        say TABOUT "$1\t0\t0\twhatever";
        next;
    }

    # Parse out the sequence name of proteins with repeats
    if(/^>(\S+)/){
        # Number of MSA encountered so far under current seqid
        $msa_count = 0; 
        $seqid = $1;
        next;
    }

    # Parser for extracting line of info for identified motif
    if(/XXX/){
        $msa_count++;
        say TABOUT "$seqid\tXXX\tXXX\t....";
        say MSAOUT ">$seqid";
    }

    # Parser for writing to MSA file
    if(/XXX/){
        say MSAOUT $1;
    }
}

exit;
