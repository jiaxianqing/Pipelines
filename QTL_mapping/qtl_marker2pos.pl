#!/usr/bin/perl -w
# jiaxianqing 2017-10-26
# alignment results to csv file.
# clustalo => alignment


use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.10.26';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($input, $posfile, $cmfile, $rate, $output, $show_help);
GetOptions(
            
            "i|input=s"            => \$input,
            "pos|pos-file=s"       => \$posfile,
            "cm|cm-file=s"         => \$cmfile,
            "rt|recombination-rate=s"       => \$rate,
            "o|output=s"           => \$output,
            "help|?"               => \$show_help,
           );

unless( !$show_help && $cmfile ) {
    print <<EOF;
$0  -- alignment results to csv file.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -i, --input <filename>
        required
    
    -pos, --pos-file <filename>
        input pos-file
    -cm, --cm-file <filename>
        input cm-file
    
    -rt, --recombination-rate [31]
        recombiantion rate cM/Mb
        
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

print STDERR ">> Start to process: $cmfile\n";

################################################################ MAIN ###################################################################
$rate = 31 if(!$rate);
my %pos = ();
open POS, "<$posfile";
while(<POS>){
    chomp;
    my @line = split /\t/, $_;
    $pos{$line[0]} = $line[2];
}
close POS;
print STDOUT "#Marker\tChrom\tPos\tLOD\n";
open CM, "<$cmfile";
while (<CM>){
    chomp;
    next if (/pos/);
    my @line = split /\t/, $_;
    if (!$pos{$line[0]}){
        my $aa = int $line[2] * 373245519 / ($rate * 100);
        $pos{$line[0]} = $aa;
    }
    print STDOUT "$line[0]\tchr$line[1]\t$pos{$line[0]}\t$line[3]\n";
}
close CM;
