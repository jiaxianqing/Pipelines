#!/usr/bin/perl -w
#
#   Author: jiaxianqing
#   Created: 2018-10-08
#   Version: 1.0.18.10.08
#   Updated: 

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.18.09.26';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my %options = ();
my %bed = ();
my ($input, $output, $show_help);
GetOptions(
            "i|input=s"              => \$input,
            
            "o|output=s"             => \$output,
            
            "help|?"               => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  transform vcf file to r/qtl format.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> -genome <genome filename>

Note:   all the range is :[a,b)

Options:
    -i, --input <filename>
        required
    
    -o, --output   <filename>
        output filename, default to STDOUT
    -?, --help
        show this help message

EOF

    exit(1);
}

$|++;


print STDERR "\n# $0 v$VERSION\n# " . (scalar localtime()) . "\n\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

print STDERR ">> Start to process: $input\n";


####################################################### Main ##############################################################
my %query = ();
open IN, "<$input";
while (<IN>){
    chomp;
    next if (/^#/);
    my @line = split /\t/, $_;
    $line[0] =~ /(.*)_(.*)/;
    my $gene = $1;
    my $type = $2;
    $query{$gene}{$type} ++;
}
close IN;

foreach my $gene (sort keys %query){
    my @out = ();
    push @out, $gene;
    foreach my $aa (sort keys %{$query{$gene}}){
        push @out, $gene."_"."$aa\t$query{$gene}{$aa}";
    }
    my $out = join "\t", @out;
    print STDOUT "$out\n";
}
