#!/usr/bin/perl -w
#   fasta2primer3.pl -- generate input format of primer3
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
my ($input, $stats, $output, $show_help);
GetOptions(
            "i|input=s"              => \$input,
            
            "s|stats=s"              => \$stats,
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
    
    -s, --stats <filename>
        output information of primer pairs

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


####################################################### prepare #############################################################

#preparation for header
open OUT, ">$stats";
print OUT "#Seq_name\tPrimer\tType\tStart\tEnd\tTm\tGC%\tAny_th\t3\'_th\tHairpin\tPrimer_seq\tProduct_len\n";
my $id = "";
my $count = 1;
my $left = 0;
open IN, "<$input";
while (<IN>){
    chomp;
    if (/PRIMER PICKING RESULTS FOR (.*)/){
        $id = $1;
        $count = 1;
    }
    my @out = ();
    push @out, $id;
    
    if($_ =~ /(LEFT PRIMER.*)/){
        my $in = $1;
        $in =~ s/PRIMER//;
        my $name = $id."_".$count."_U";
        push @out, $name;
        my @line = split /\s+/, $in;
        foreach my $l (@line){
            push @out, $l;
        }
        my $out = join "\t", @out;
        print OUT "$out\n";
        print STDOUT ">$name\n$line[-1]\n";
        $left = $line[1];
    }elsif($_ =~ /(RIGHT PRIMER.*)/){
        my $in = $1;
        $in =~ s/PRIMER//;
        my $name = $id."_".$count."_L";
        push @out, $name;
        my @line = split /\s+/, $in;
        foreach my $l (@line){
            push @out, $l;
        }
        push @out, $line[1] - $left + 1;
        my $out = join "\t", @out;
        print OUT "$out\n";
        print STDOUT ">$name\n$line[-1]\n";
        $count ++;
    }
}
