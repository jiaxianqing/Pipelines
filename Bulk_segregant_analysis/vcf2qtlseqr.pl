#!/usr/bin/perl -w
# jiaxianqing 2019-01-07
# convert VCF file format to QTLseqr input format
# only keep thess fields in VCF file: CHROM, POS, REF, ALT, AD, DP, GQ
# https://github.com/bmansfeld/QTLseqr

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.01.07';
my $HEADER  = "Version: $VERSION $CMDLINE\n";

my %options = ();
my %bed = ();
my ($input, $groupfile, $method, $maf, $mdp, $output, $show_help);
GetOptions(
            "i|input=s"                          => \$input,
            "g|group=s"                          => \$groupfile,
            "m|merge-method=s"                   => \$method,
            "maf|minor-allele-frequency=s"       => \$maf,
            "mdp|minor-total-depth=s"            => \$mdp,
            "o|output=s"                         => \$output,
            "help|?"                             => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  -- alignment results to csv file.

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    -i, --input <filename>
        required
    
    -g, --group <filename>
        group file name, only keep samples listed in this file.
        e.g.
        #Sample Group
        DE1 Early
        DN2 Late
    
    -m, --merge-method <AF or DP>
        how to merge samples together, default: DP:
        GT:AD:DP 0/1:21,23:44 0/1:21,23:44 0/1:21,23:44 0/1:21,23:44 
        1) AF: merge by allele frequency, 0/1:4,4:8
        2) DP: merge by depth, 0/1:84,92:176
    
    -maf, --minor-allele-frequency <number>
        you can assign a filter by MAF. [0]
    -mdp, --minor-total-depth <number>
        you can assign a filter by MDP. [0]
    
    -o, --output   <filename>
        output filename, default to STDOUT
        CHROM	POS	REF	ALT	93-11.AD	93-11.DP	93-11.GQ	DE10.AD	DE10.DP	DE10.GQ
        chr01	1023	TA	T	13,5	18	99	8,4	12	99
        chr01	1044	TA	T	12,9	21	99	9,7	16	99
        chr01	1058	T	TA	13,12	25	99	10,10	20	99
        
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
#CHROM	POS	REF	ALT	93-11.AD	93-11.DP	93-11.GQ	93-11.PL	DE10.AD	DE10.DP	DE10.GQ	DE10.PL
#chr01	1023	TA	T	13,5	18	99	149,0,429	8,4	12	99	127,0,256
#CHROM POS REF ALT AD_REF.SRR834927 AD_ALT.SRR834927 GQ.SRR834927 AD_REF.SRR834931 AD_ALT.SRR834931 GQ.SRR834931
if(!$method){
    $method = "DP";
}


my %sample = ();
my %group = ();
open G, "<$groupfile";
while (<G>){
    chomp;
    next if (/^#/);
    my @line = split /\t/, $_;
    $sample{$line[1]} ++;
    $group{$line[1]}{$line[0]} = 1;
}
close G;


my %pos = ();
my %num = ();
open IN, "<$input";
while(<IN>){
    chomp;
    next if(/^##/);
    my @line = split /\t/, $_;
    if(/^#C/){
        for (my $i = 9; $i <= $#line; $i ++){
            $pos{$line[$i]} = $i;
            $num{$i} = $line[$i];
        }
        my @out = ();
        foreach my $g (sort keys %group){
            push @out, "AD_REF.$g\tAD_ALT.$g\tGQ.$g";
        }
        my $out = join "\t", @out;
        print STDOUT "CHROM\tPOS\tREF\tALT\t$out\n";
    }else{
        next if ($line[4] =~ /,/);
        my @out = ();
        push @out, $line[0];
        push @out, $line[1];
        push @out, $line[3];
        push @out, $line[4];
        #GT:AD:DP:GQ:PL 0/0:32,1:33:84:0,84,1100
        #AD, DP, GQ
        my $flag = 0;
        if($method eq "AF"){
            foreach my $g (sort keys %group){
                my @ad = ();
                $ad[0] = 0;
                $ad[1] = 0;
                my $gq = 0;
                my $s = 0;
                for (my $i = 9; $i <= $#line; $i++){
                    next if ($line[$i] =~ /\.\/\./);
                    if ($group{$g}{$num{$i}}){
                        my @info = split /:/, $line[$i];
                        my @gt = split/\//, $info[0];
                        $ad[$gt[0]] ++;
                        $ad[$gt[1]] ++;
                        $gq += $info[3];
                        $s ++;
                    }
                }
                if ($s < $mdp){
                    $flag = 1;
                    next;
                }
                my $total_allele = $ad[0] + $ad[1];
                if ($ad[0]/$total_allele < $maf or $ad[1]/$total_allele < $maf){
                    $flag = 1;
                    next;
                }
                my $avr_gq = sprintf('%.f', $gq / $s);
                push @out, "$ad[0]\t$ad[1]";
                #push @out, $total_allele;
                push @out, $avr_gq;
            }
            
        }elsif($method eq "DP"){
            foreach my $g (sort keys %group){
                my @ad_sum = ();
                $ad_sum[0] = 0;
                $ad_sum[1] = 0;
                my $gq = 0;
                my $dp = 0;
                my $s = 0;
                for (my $i = 9; $i <= $#line; $i++){
                    next if ($line[$i] =~ /\.\/\./);
                    if ($group{$g}{$num{$i}}){
                        my @info = split /:/, $line[$i];
                        my @ad = split/,/, $info[1];
                        $ad_sum[0] += $ad[0];
                        $ad_sum[1] += $ad[1];
                        $dp += $info[2];
                        $gq += $info[3];
                        $s ++;
                    }
                }
                #print STDOUT "$dp\n";
                if ($dp < $mdp){
                    $flag = 1;
                    next;
                }
                if ($ad_sum[0]/$dp < $maf or $ad_sum[1]/$dp < $maf){
                    $flag = 1;
                    next;
                }
                my $avr_gq = sprintf('%.f', $gq / $s);
                push @out, "$ad_sum[0]\t$ad_sum[1]";
                push @out, $avr_gq;
            }
        }
        next if ($flag == 1);
        my $out = join "\t", @out;
        print STDOUT "$out\n";
    }
}
close IN;