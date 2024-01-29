#!/usr/bin/perl -w

# Xianqing Jia
# 2019-03-30
# trans gtf format to table

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.03.30';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my ($input, $list, $output, $sample, $show_help);
GetOptions(
            "i|input=s"                => \$input,
            "l|list=s"                 => \$list,
            "s|samples=s"              => \$sample,
            "o|output=s"               => \$output,
            "help|?"                   => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  -- vcf2matrix and regenotype.
    PS. remove indels

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    [file options]
    -i, --input <filename>
        required. vcf file or table file(from vcf2table.pl)
    -l, --list <filename>
        query list
    -s, --sample <name>
        
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

################################################################ MAIN ###################################################################
my %sample = ();
my %tpm = ();
open IN, "<$input";
while (<IN>){
    chomp;
    open SAM, "<$_";
    while (<SAM>){
        chomp;
        next if (/^#/);
        my @line = split /\t/, $_;
        $tpm{$line[2]}{$line[0]} = $_;
        $sample{$line[0]} = 1;
    }
    close SAM;
}
close IN;

my @sam = ();
foreach my $s (sort keys %sample){
    push @sam, $s;
}
my $sam_out = join "\t", @sam;
print STDOUT "#Gene\tTranscript\t$sam_out\n";


foreach my $trans (sort keys %tpm){
    my @aa = split /\./, $trans;
    my @out = ();
    my $gene = "";
    foreach my $s (sort keys %sample) {
        my @bb = split /\t/, $tpm{$trans}{$s};
        $gene = $bb[1];
        push @out, $bb[-1];
    }
    my $out = join "\t", @out;
    print STDOUT "$gene\t$trans\t$out\n";
}

