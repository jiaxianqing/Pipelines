#!/usr/bin/perl -w
#
#   find_low_copy_genes.pl -- find low-copy (1-2) gene orthologs
#
#   Author: jia xianqing
#   Created: 2017-06-11
#   Updated: 
#   Version: 1.0.17.06.11
#
#   Change logs:
#   Version 17.06.11 2017-06-11: The initial version.

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.06.11';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my %options = ();
my %bed = ();
my ($seqdir, $suffix, $lowcopygroup, $outdir, $show_help);
GetOptions(
            "sd|seqdir=s"              => \$seqdir,
            "sf|suffix=s"              => \$suffix,
            "lcg|lowcopygroup=s"       => \$lowcopygroup,
            
            "od|outdir=s"              => \$outdir,           
            "help|?"                   => \$show_help,
           );

unless( !$show_help && $seqdir ) {
    print <<EOF;
$0  -- count site types in certain region.
    PS. surpport multi-type

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -sd, --seqdir <dirname>
        required
    -sf, --suffix <string>
        required, e.g. fa, fas, pep
    
    -lcg, --lowcopygroup   <filename>
        output ortholog group details
    -od, --oudtir   <dirname>
        output directory
    
    -?, --help
        show this help message

EOF

    exit(1);
}

$|++;


print STDERR "\n# $0 v$VERSION\n# " . (scalar localtime()) . "\n\n";


print STDERR ">> Start to process: $seqdir\n\n";
print STDERR ">> Wraning: missing sequences will be shown flow: \n\n";
################################################################ MAIN ###################################################################

my %seq = ();
my @file = glob "$seqdir/*.$suffix";
foreach my $file (@file){
    open IN, "<$file";
    $file =~ /$seqdir\/(.*?)\.$suffix/;
    my $name = $1;
    my $flag = "";
    while (<IN>){
        chomp;
        if(/^>(.*?)\s+/){
            $flag = $1;
            next;
        }elsif(/^>(.*)/){
            $flag = $1;
        }else{
            $seq{$name}{$flag} .= $_;
        }
    }
    close IN;
}

my %order = ();
open IN, "<$lowcopygroup";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    if(/^Ortho_group/){
        for (my $i = 1; $i <= $#line; $i ++){
            $order{$i} = $line[$i];
        }
        next;
    }else{
        open OUT, ">$outdir/$line[0].fa";
        for (my $i = 1; $i <= $#line; $i ++){
            if($seq{$order{$i}}{$line[$i]}){
                print OUT ">$line[$i]\n$seq{$order{$i}}{$line[$i]}\n";
            }else{
                print ">$line[$i]\n----\n"
            }
        }
        close OUT;
    }
}
close IN;


