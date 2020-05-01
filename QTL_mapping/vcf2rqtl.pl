#!/usr/bin/perl -w



#   drop_site_into_blocks.pl -- find sites which are droped into certin blocks
#
#
#   Author: jiaxianqing
#   Created: 2016-06-22
#   Version: 1.0.16.06.22
#   Updated: 1.1.16.07.05 change calculated method of genetic map distance(cM):
#            rice genome total length: 373245519
#            average CO number: 31, 3100cM
#            patchwork: 3100/373245519 cM/bp
#            1.2.18.01.25 add output option: markers position
#            1.2.18.03.14 change extracted method of phenotype name(get from title of phenotype file)
#                         add recombination times (-rt, --rc-times) options
#            1.3.18.05.10 add symbol (-sy, --symbol) option for SC tag in the vcf file
#            1.4.18.07.25 add option -- file of removed samples.
use strict;
use Getopt::Long;



my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.4.18.07.25';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my %options = ();
my %bed = ();
my ($vcf, $genome, $pheno, $remove, $rctimes, $removefile, $output, $marker, $symbol, $show_help);
GetOptions(
            "v|vcf=s"              => \$vcf,
            
            "p|phenotype=s"        => \$pheno,
            "g|genome=s"           => \$genome,
            "r|remove=s"           => \$remove,
            "rf|remove-file=s"     => \$removefile,
            "rt|rc-times=i"        => \$rctimes,
            "sy|symbol=s"          => \$symbol,
            "o|output=s"             => \$output,
            "m|maker=s"            => \$marker,
            "help|?"               => \$show_help,
           );

unless( !$show_help && $vcf ) {
    print <<EOF;
$0  transform vcf file to r/qtl format.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> --phenotype <phenotype filename> --genome <genome filename>

Note:   all the range is :[a,b)

Options:
    -v, --vcf <filename>
        required
    -p, --phenotype <filename>
        required
    
    -g, --genome <filename>
        required
    -rc, --rc-times <int>
        recombination times in the whole genome. [34]
    
    -r, --remove <sample name>
        separated by comma
    -rf, --remove-file <file name>
        list removed sample names in a file, one name per line.
        
    -sy, --symbol <symbol name>
        symbol of parent genotypes, SC tag in the vcf, separated by comma. [A,B]
    -m, --marker <filename>
        output markers position and corresponding table
    
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

print STDERR ">> Start to process: $vcf\n\n";


####################################################### prepare #############################################################

##prepare for transform

## read genome file to get chromesome length
my %genome = ();
my $gsize = 0;
open IN, "<$genome";
while (<IN>){
    chomp;
    my @line = split /\t/, $_;
    $genome{$line[0]} = $line[1];
    $gsize += $line[1];
}
close IN;

## read phenotype file
my %pheno = ();
my @ph = ();
open IN, "<$pheno";
while (<IN>){
    chomp;
    my @line = split /\s+/, $_;
    if(/^#/){
        push @ph, $line[1];
        next;
    }
    $pheno{$line[0]} = $line[1];
    
}
close IN;
push @ph, "-";
push @ph, "-";
## get removed sample;
my %remove = ();
if($remove){
    my @re = split /,/, $remove;
    foreach my $r (@re){
        $remove{$r} = 1;
    }
}
if($removefile){
    open RF, "<$removefile";
    while (<RF>){
        chomp;
        $remove{$_} = 1;
    }
    close RF;
}
## get recombination times
if(!$rctimes){
    $rctimes = 34;
}
## get symbol of parent genotypes
my $sya = "A";
my $syb = "B";

if($symbol){
    my @sy = split /,/, $symbol;
    $sya = $sy[0];
    $syb = $sy[1];
}
####################################################### MAIN #############################################################
open M, ">$marker";
## read vcf file
my $chr = "";
my $count = 1;

my %sample = ();
open IN, "<$vcf";
while (<IN>) {
    chomp;
    next if (/^##/);
    s/$sya\/$sya/AA/g;
    s/$sya\/$syb/AB/g;
    s/$syb\/$syb/BB/g;
    s/.\/./NA/g;
    
    my @out = ();
    my @line = split /\t/, $_;
    if (/^#C/){
        for (my $i = 9; $i <= $#line; $i ++){
            $sample{$i} = $line[$i];
            if ($remove{$line[$i]}){
                next;
            }else{
                push @ph, $pheno{$line[$i]};
            }
            
        }
        my $out = join ",", @ph;
        print STDOUT "$out\n";
    }else{
        if ($chr eq $line[0]){
            $count ++;
        }else{
            $count = 1;
            $chr = $line[0];
        }
        my $a = "$line[0]" . "S" . "$count";
        push @out, $a;
        push @out, $line[0];
        my $b = $line[1] * 100 * $rctimes / $gsize;
        push @out, $b;
        
        for (my $i = 9; $i <= $#line; $i ++){
            if ($remove{$sample{$i}}){
                next;
            }else{
                push @out, $line[$i];
            }
        }
        my $out = join ",", @out;
        print STDOUT "$out\n";
        print M "$a\t$b\t$line[1]\n";
    }
}
close IN;
close M;