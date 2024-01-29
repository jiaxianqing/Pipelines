#!/usr/bin/perl -w

# Xianqing Jia
# 2019-03-30

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.03.30';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my ($input, $output, $show_help);
GetOptions(
            "i|input=s"                => \$input,
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

my $lc = 0;
my %all = ();
open IN, "<$input";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    if(/^#/){
        s/Transcript\t//;
        print STDOUT "$_\n";
        $lc = $#line;
        next;
    }else{
        for (my $i = 2; $i <= $#line; $i ++){
            $all{$line[0]}{$i} += $line[$i];
        }
    }
}
close IN;

foreach my $g (sort keys %all){
    my @out = ();
    for (my $s = 2; $s <= $lc; $s ++) {
        push @out, $all{$g}{$s};
    }
    my $out = join "\t", @out;
    print STDOUT "$g\t$out\n";
}
