#!/usr/bin/perl -w
#
#   count_markers.pl -- count marker number in certain region
#
#   Author: jia xianqing
#   Created: 2017-06-11
#   Updated: 
#   Version: 1.0.17.06.11
#
#   Change logs:
#   Version 1.0.0 15/11/15: The initial version.
use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.06.11';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my %options = ();
my %bed = ();
my ($indir, $suffix, $output, $lowcopygroup, $show_help);
GetOptions(
            
            "id|indir=s"           => \$indir,
            "sf|suffix=s"          => \$suffix,
            "lcg|lowcopygroup=s"   => \$lowcopygroup,
            
            "o|output=s"            => \$output,
            "help|?"               => \$show_help,
           );

unless( !$show_help && $indir ) {
    print <<EOF;
$0  -- count site types in certain region.
    PS. surpport multi-type

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -i, --input <filename>
        required
    -g, --genome <filename>
        genome length file(.fai)
    -t, times <int>
        simulate times
    
    -o, --output   <filename>
        output filename, default to STDOUT
    
    -?, --help
        show this help message

EOF

    exit(1);
}

$|++;


print STDERR "\n# $0 v$VERSION\n# " . (scalar localtime()) . "\n\n";


print STDERR ">> Start to process: $indir\n\n";

################################################################ MAIN ###################################################################

my %seq = ();
my @file = glob "$indir/*.$suffix";
foreach my $file (@file){
    open IN, "<$file";
    my $flag = 0;
    while (<IN>){
        chomp;
        if(/^>/){
            $flag ++;
            #print "$flag\n";
        }else{
            $seq{$flag} .= $_;
        }
    }
    close IN;
}

my %order = ();
open IN, "<$lowcopygroup";
while (<IN>){
    chomp;
    my @line = split /\t/;
    if(/^Ortho_group/){
        s/.pep//;
        for (my $i = 1; $i <= $#line; $i ++){
            $order{$i} = $line[$i];
        }
        next;
    }
}
open STDOUT, ">$output";
foreach my $s (sort keys %seq){
    print STDOUT ">$order{$s}\n$seq{$s}\n";
}
close STDOUT;