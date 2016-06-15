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
my $seqid_count;
while(<INPUT>){
    if(/^>(\S+)/){
        $seqid_count = 0; 
        $seqid = $1;
        next;
    }

    if(/not found.*>(\S+)/) {
        say TABOUT "$1\t0\t0\twhatever";
        next;
    }

    # other parsers
}

exit;
