#!/usr/bin/perl -w

# Xianqing Jia
# 2019-01-10


use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.10.16';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($input, $site_file, $difference, $output, $show_help);
GetOptions(
            "i|input=s"            => \$input,
            "s|site=s"           => \$site_file,
            "d|difference=s"        => \$difference,
            "o|output=s"         => \$output,
            "help|?"             => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0 -- extend High-LOD sites to intervals.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -i, --input <filename>
        required
    -s, --site <filename>
        informations of all sites.
        #Marker	Chrom	Pos	LOD
        chr01S1	chr01	1117	1.49923383040385
        chr01S2	chr01	1131	1.49708026125734
        chr01S3	chr01	1132	0.830738711658086
    
    -d, --difference [1.8]
        changeable intervals by 1.8 region.
        for a site, LOD equal n, region: n-1.8 (left), n, n-1.8(right)
        
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
my %site = ();
my %marker = ();
my %pos = ();
my %order = ();
my %lod = ();
my %chr = ();
my $count = 1;
open IN, "<$site_file";
while (<IN>){
    chomp;
    next if (/^#/);
    my @line = split /\t/, $_;
    $site{$line[0]} = $line[2];
    $chr{$line[0]} = $line[2];
    $lod{$line[0]} = $line[3];
    $order{$count} = $line[0];
    $pos{$line[0]} = $count;
    $count ++;
}
close IN;

print STDOUT "#Marker\tChrom\tPos\tLOD\tStart\tEnd\tLength\tAver_LOD\n";
open IN, "<$input";
while (<IN>){
    chomp;
    next if (/^#/);
    my @line = split /\t/, $_;
    my @out = ();
    push @out, $_;
    my $poshigh = $pos{$line[0]};
    my $leftfinal = $line[0];
    my $mc = 1;
    my $sum = $lod{$line[0]};
    for (my $i = 1; $i <= $count; $i ++){
        my $leftpos = $poshigh - $i;
        my $leftmar = $order{$leftpos};
        my @l = split /S/, $leftmar;
        if($l[0] ne $line[1]){
            last;
        }
        if($lod{$leftmar} >= $line[3] - $difference){
            $leftfinal = $leftmar;
            $mc ++;
            $sum += $lod{$leftfinal};
        }else{
            last;
        }
    }
    my $rightfinal = $line[0];
    for (my $i = 1; $i <= $count; $i ++){
        my $rightpos = $poshigh + $i;
        my $rightmar = $order{$rightpos};
        my @r = split /S/, $rightmar;
        if($r[0] ne $line[1]){
            last;
        }
        if($lod{$rightmar} >= $line[3] - $difference){
            $rightfinal = $rightmar;
            $mc ++;
            $sum += $lod{$rightfinal};
        }else{
            last;
        }
    }
    my $averlod = sprintf "%.2f", $sum / $mc;
    my $length = $site{$rightfinal} - $site{$leftfinal} + 1;
    push @out, "$site{$leftfinal}\t$site{$rightfinal}\t$length\t$averlod";
    my $out = join "\t", @out;
    print STDOUT "$out\n";
}