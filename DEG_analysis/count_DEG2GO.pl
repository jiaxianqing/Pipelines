#!/usr/bin/perl -w
# jiaxianqing 2017-03-14

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.10.16';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($input, $uplist, $downlist, $output, $show_help);
GetOptions(
            
            "i|input=s"          => \$input,
            "u|uplist=s"         => \$uplist,
            "d|downlist=s"       => \$downlist,
            "o|output=s"         => \$output,
            "help|?"             => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  -- vcf2matrix and regenotype.
    PS. remove indels

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:

    
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
my %up = ();
open IN, "<$uplist";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    $up{$line[0]} = 1;
}
close IN;

my %down = ();
open IN, "<$downlist";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    $down{$line[0]} = 1;
}
close IN;

open GO, "<$input";
while (<GO>){
    chomp;
    if (/^GO_acc/){
        print STDOUT "$_\tUp_number\tUp_gene\tDown_number\tDown_gene\n";
        next;
    }else{
        my @line = split /\t/, $_;
        $line[-1] =~ s/^\/\///;
        $line[-1] =~ s/ //g;
        my @gene = split /\/\//, $line[-1];
        my $up = 0;
        my @upgene = ();
        my $down = 0;
        my @downgene = ();
        foreach my $g (@gene){
            if ($up{$g}){
                push @upgene, $g;
                $up ++;
            }elsif($down{$g}){
                push @downgene, $g;
                $down ++;
            }
        }
        my $upout = join ";", @upgene;
        my $downout = join ";", @downgene;
        print STDOUT "$_\t$up\t$upout\t$down\t$downout\n";
    }
}
close GO;
