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
my ($lowcopylist, $orthofile, $genecount, $lowcopygroup, $output, $show_help);
GetOptions(
            "gc|genecount=s"     => \$genecount,
            "og|orthogroups=s"   => \$orthofile,
            
            "lcl|lowcopylist=s"       => \$lowcopylist,
            "lcg|lowcopygroup=s"       => \$lowcopygroup,
            
            "help|?"               => \$show_help,
           );

unless( !$show_help && $genecount ) {
    print <<EOF;
$0  -- find low-copy (1-2) gene orthologs.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -gc, --genecount <filename>
        required
    -og, --orthogroups <filename>
        required
        
    -lcl, --lowcopylist <filename>
        output list of ortholog group names with low copy genes
    
    -lcg, --lowcopygroup   <filename>
        output ortholog group details
    
    -?, --help
        show this help message

EOF

    exit(1);
}

$|++;


print STDERR "\n# $0 v$VERSION\n# " . (scalar localtime()) . "\n\n";


print STDERR ">> Start to process: $genecount\n\n";

################################################################ MAIN ###################################################################

my %list = ();
open IN, "<$genecount";
open OUT, ">$lowcopylist";
while (<IN>){
    chomp;
    #next if (/Amborella/);
    next if (/^\t/);
    my @line = split /\t/, $_;
    my $flag = 0;
    for (my $i = 1; $i < $#line; $i ++){
        if ($line[$i] eq 0 or $line[$i] > 2){
            $flag = 1;
        }
    }
    if ($flag == 0){
        print OUT "$line[0]\n";
        $list{$line[0]} = 1;
        #print OUT "$_\n";
    }
}
close IN;
close OUT;


open IN, "<$orthofile";
open OUT, ">$lowcopygroup";
while (<IN>){
    chomp;
    if (/^\t/){
        my $out = "Ortho_group" . $_;
        print OUT "$out\n";
    }
    
    my @line = split /\t/, $_;
    if ($list{$line[0]}){
        my @out = ();
        push @out, $line[0];
        for (my $i = 1; $i <= $#line; $i ++){
            if ($line[$i] =~ /, /){
                my @aa = split /, /, $line[$i];
                push @out, $aa[0];
            }else{
                push @out, $line[$i];
            }
        }
        my $out = join "\t", @out;
        print OUT "$out\n";
    }
}
close IN;
close OUT;
