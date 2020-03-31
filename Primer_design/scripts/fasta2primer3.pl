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
my ($input, $genome, $task, $number, $target_region, $size_region, $output, $show_help);
GetOptions(
            "i|input=s"                       => \$input,
            
            "t|task=s"                        => \$task,
            "tr|target_region=s"              => \$target_region,
            "n|number=s"                      => \$number,
            "s|size-region=s"                 => \$size_region,
            
            "g|genome=s"           => \$genome,
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
    
    -t, --task <name>
        required, OPTIONS: generic, pick_primer_list, check_primers, pick_pcr_primers et.al.
    -tr, --target_region <start,length>
        required
    -s, --size-region <min-max>
        required
    -n, --number <num>
        required, numbers of candidated primer pair
    
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


####################################################### MAIN #############################################################
#preparation for header

print STDOUT "PRIMER_TASK=$task\n";
print STDOUT "PRIMER_PRODUCT_SIZE_RANGE=$size_region\n";
print STDOUT "SEQUENCE_TARGET=$target_region\n";
print STDOUT "PRIMER_NUM_RETURN=$number\n=\n";


my %id = ();
open IN, "<$input";
while (<IN>){
    if(/^>(.*)/){
        my $name = $1;
        chomp(my $seq = <IN>);
        my $prefix = substr $seq, 0, 5;
        if($id{$name}){
            $name = $name . "\.". $prefix;
        }else{
            $id{$name} = 1;
        }
        print STDOUT "SEQUENCE_ID=$name\nSEQUENCE_TEMPLATE=$seq\n=\n";
    }
}
close IN;

