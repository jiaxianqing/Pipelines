#!/usr/bin/perl -w
# jiaxianqing 2018-09-09

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.10.26';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($input, $name, $output, $symbol, $show_help);
GetOptions(
            
            "i|input=s"            => \$input,
            "n|name=s"             => \$name,
            
            "o|output=s"            => \$output,
            "s|symbol=s"               => \$symbol,
            "help|?"               => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  -- alignment results to csv file.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -i, --input <filename>
        required
    
    -d, --domain
        domain file.
    
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
my %list = ();
open IN, "<$symbol";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    $list{$line[0]} = $line[1];
}
close IN;

open IN, "<$input";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    #$line[1] =~ s/\.v5.5//g;
    if ($list{$line[1]}){
        print STDOUT "$name\t$list{$line[1]}\t$_\n";
    }else{
        print STDOUT "$name\tNA\t$_\n";
    }
}
close IN;

