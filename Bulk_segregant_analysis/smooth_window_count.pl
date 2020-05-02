#!/usr/bin/perl -w
# jiaxianqing 2017-03-14
# for re-genotyping the vcf file


use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.17.10.16';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($bed_file, $site_file, $type, $method, $site_info, $output, $show_help);
GetOptions(
            "b|bed=s"            => \$bed_file,
            "s|site=s"           => \$site_file,
            "m|methods=s"        => \$method,
            "o|output=s"         => \$output,
            "help|?"             => \$show_help,
           );

unless( !$show_help && $bed_file ) {
    print <<EOF;
$0  -- vcf2matrix and regenotype.
    PS. remove indels

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -b, --bed <filename>
        required
    -s, --site <filename>
        count site number of each block and stats for markers by "-m, --methods".
        #chr    pos value1  value2
        chr01   1   5   0.1
        chr01   4   6   0.2
    -m, --methods [average]
        average: average value of all sites.
        proportion: length proportion of the whole window.
        
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

print STDERR ">> Start to process: $bed_file\n";

################################################################ MAIN ###################################################################
#my $TIME = scalar localtime();
#print STDOUT "##$TIME $HEADER";

$method = "average" if (!$method);

my %type = ();
my %sites = ();
open SITE, "<$site_file";
while (<SITE>){
    chomp;
    next if (/^#/);
    my @line = split /\s+/, $_;
    $sites{$line[0]}{$line[1]} = $line[2];
}
close SITE;

print STDOUT "Chrom\tStart\tEnd\tLength\tBin\tSite_count\tStats\n";
my $flag = 1;
my $chr = "NA";
open BED, "<$bed_file";
while (<BED>){
    chomp;
    next if (/^#/);
    my $in = $_;
    my @line = split /\t/, $in;
    my @output = ();
    my $length = $line[2] - $line[1];
    push @output, $in;
    push @output, $length;
    if($chr ne $line[0]){
        $flag = 1;
        $chr = $line[0];
    }else{
        $flag ++;
    }
    push @output, $flag;
    if($method eq "average"){
        my @stats = &AverageStats ($in, %sites);
        #my $aa = &Count ($in, %sites);
        #push @output, "$aa";
        push @output, "$stats[0]\t$stats[1]";
    }
    my $out = join "\t", @output;
    print STDOUT "$out\n";
}
close BED;
################################################################ SUB ###################################################################
sub Count {
    my ($block, %site) = @_;
    my @block = split /\s+/, $block;
    my @num = grep {$_ > $block[1] and $_ <= $block[2]} sort { $a <=> $b } keys %{$site{$block[0]}};
    #@num = sort @num;
    my $aa = join ",", @num;
    return $aa;
    
}

sub AverageStats {
    my ($block, %site) = @_;
    my @block = split /\s+/, $block;
    my @num = grep {$_ > $block[1] and $_ <= $block[2]} sort keys %{$site{$block[0]}};
    my $sum = 0;
    my @out = ();
    if(@num){
        foreach my $n (@num){
            $sum += $site{$block[0]}{$n};
        }
        my $count = $#num + 1;
        my $average = $sum / $count;
        push @out, $count;
        push @out, $average;
    }else{
        push @out, "0";
        push @out, "0";
    }
    return @out;
}
