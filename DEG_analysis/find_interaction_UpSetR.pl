#!/usr/bin/perl -w

# Xianqing Jia
# 2019-03-30

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.03.30';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my ($indir, $output, $fc, $suffix, $show_help);
GetOptions(
            "i|indir=s"                => \$indir,
            "s|suffix=s"               => \$suffix,
            "log2fc|log2fc=s"          => \$fc, 
            "o|output=s"               => \$output,
            "help|?"                   => \$show_help,
           );

unless( !$show_help && $indir ) {
    print <<EOF;
$0  -- vcf2matrix and regenotype.
    PS. remove indels

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    [file options]
    -i, --indir <filename>
        required. vcf file or table file(from vcf2table.pl)

    -s, --suffix <name>
        
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
print STDERR ">> Start to process: $indir\n";

################################################################ MAIN ###################################################################

my %all = ();
my %all_gene = ();
my @file = glob "$indir/*.$suffix";
foreach my $file (@file){
    $file =~ /.*\/(.*?).DEseq2/;
    my $samp = $1;
    print STDERR "$samp\n";
    open IN, "<$file";
    while (<IN>){
        chomp;
        next if (/^Row.names/);
        my @line = split /\t/, $_;
        next if (abs($line[2]) <= $fc);
        $all{$samp}{$line[0]} = 1;
        $all_gene{$line[0]} = 1;
    }
    close IN;
}

my @out = "";
foreach my $file (@file){
    #print OUT "$file\n";
    $file =~ /.*\/(.*?).DEseq2/;
    my $samp = $1;
    push @out, $samp;
}
my $op = join "\t", @out;
print STDOUT "Gene"."$op\n";


foreach my $g (sort keys %all_gene){
    my @out1 = ();
    foreach my $t (@out){
        if ($all{$t}{$g}) {
            push @out1, "1";
        }else{
            push @out1, "0";
        }
    }
    my $op = join "\t", @out1[1..$#out1];
    print STDOUT "$g\t$op\n";

}
