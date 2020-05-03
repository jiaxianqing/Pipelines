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
my ($input, $threshold, $output, $show_help);
GetOptions(
            "i|input=s"                => \$input,
            "t|threshold=s"            => \$threshold,
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
    -t, --threshold <filename>
        threshold value
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
print STDOUT "#Gene\tTranscript\tChrome\tStart\tEnd\tAvg_cov\tFPKM\tTPM\n";
open IN, "<$input";
while (<IN>){
    chomp;
    next if (/^#/);
    my @line = split /\t/, $_;
    if ($line[2] eq "transcript"){
        $line[-1] =~ /gene_id "(.*)"; transcript_id/;
        $line[-1] =~ /gene_id "(.*)"; transcript_id "(.*)";.*; cov "(.*)"; FPKM "(.*)"; TPM "(.*)";/;
        #print STDOUT "$1\t$line[0]\t$line[3]\t$line[4]\n";
        print STDOUT "$1\t$2\t$line[0]\t$line[3]\t$line[4]\t$3\t$4\t$5\n";
    }
}
close IN;
