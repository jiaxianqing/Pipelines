#!/usr/bin/perl -w

# Xianqing Jia
# 2019-11-26
# assign experimental design to each dir

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.11.26';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my ($input, $count_dir, $output_dir, $show_help);
GetOptions(
            "i|input=s"                => \$input,
            "cd|count_dir=s"           => \$count_dir,
            "od|output_dir=s"          => \$output_dir,
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
print STDERR ">> Start to process: $input\n";

################################################################ MAIN ###################################################################
open IN, "<$input";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    if (/^Group/){
        mkdir "$output_dir/$line[0]";
        open OUT1, ">$output_dir/$line[0]/phenotype.tab";
        print OUT1 "sampleID\tgroup\n";
        close OUT1;
        open OUT2, ">$output_dir/$line[0]/samples.list";
        close OUT2;
    }
}
close IN;

my $flag = "";
open IN, "<$input";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    if (/^Group/){
        mkdir "$output_dir/$line[0]";
        $flag = $line[0];
    }else{
        open OUT1, ">>$output_dir/$flag/phenotype.tab";
        print OUT1 "$_\n";
        close OUT1;
        open OUT2, ">>$output_dir/$flag/samples.list";
        print OUT2 "$count_dir/$line[0]\n";
        close OUT2;
    }
}
close IN;


