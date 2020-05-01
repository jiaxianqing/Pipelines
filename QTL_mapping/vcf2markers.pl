#!/usr/bin/perl -w

# Xianqing Jia
# 2019-01-02
# find markers: 1)parents-only, only find different sites in parental genotypes;
#               2)single-parent, differences occured between one parent and progeny;
#               3)no-parents, differences only come from progeny;

# Updated: 2019-01-10, added -mmp, --max-missing-proportion
# Updated: 2019-01-11, changed allele method, now contain parents

use strict;
use Getopt::Long;

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.19.01.11';
my $HEADER  = "Version: $VERSION $CMDLINE";

my %options = ();
my %bed = ();
my ($input, $segregat, $outsymbol, $maf, $mmp, $method, $removesamp, $removefile, $maxallele, $minallele, $output, $show_help);
GetOptions(
            "i|input=s"            => \$input,
            
            "s|segregating=s"        => \$segregat,
            "os|out-symbols=s"     => \$outsymbol,
            "m|method=s"           => \$method,
            
            "r|remove-sample=s"    => \$removesamp,
            "rf|remove-sample-file=s"    => \$removefile,
            
            "maxa|max-alleles=s"    => \$maxallele,
            "mina|min-alleles=s"    => \$minallele,
            
            "maf|minor-allele-frequence=s"    => \$maf,
            "mmp|max-missing-proportion=s"  => \$mmp,
            
            "o|output=s"           => \$output,
            "help|?"               => \$show_help,
           );

unless( !$show_help && $input ) {
    print <<EOF;
$0  -- vcf2matrix and regenotype.
    PS. remove indels

Version: $VERSION

Usage:   perl $0 --vcf <vcf filename> [options]

Options:
    [file options]
    -i, --input <filename>
        required. vcf file or table file(from vcf2table.pl)
        
    -o, --output   <filename>
        output filename, default to STDOUT
    
    -s, --segregating <parent1,parent2>
        assign which are parents, separated by comma
    -os, --out-symbols <symbol1,symbol2>
        symbols for output, corresponding with segregating, separated by comma.
    
    -m, --method <method name>
        different method to find markers. [m1]
        m1, parents-only, only find different sites in parental genotypes;
        m2, single-parent, differences occured between one parent and progeny;
        m3, no-parents, differences only come from progeny;
        #site numbers: m3 > m2 > m1
    
    -r, --remove-sample <sample name>
        removed samples, separated by comma
    -rf, --remove-sample-file <filename>
        list of removed samples, one sample each line.
    
    -maxa, --max-alleles <number>
        max-alleles number. [2]
    -mina, --min-alleles <number>
        min-alleles number. [2]
    
    -maf, --minor-allele-frequence <number>
        it is required if "method" option is assigned as m2, m3 or all.
    -mmp, --max-missing-proportion <number>
        max number of samples proportion which are with missing genotype, <= mmp. [0.5]
        
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

# define parents lines

$mmp = 0.5 if (!$mmp);
$method = "m1" if (!$method);
$minallele = 2 if (!$minallele);
$maxallele = 2 if (!$maxallele);

my %symbol = ();
my @p = split /,/, $segregat;
my @s = split /,/, $outsymbol;
$symbol{$p[0]} = $s[0];
$symbol{$p[1]} = $s[1];
# add annotation in header of VCF file.
my $a1 = "##INFO=<ID=$s[0],Number=1,Type=String,Description=\"Corresponding variant of the symbol\">";
my $a2 = "##INFO=<ID=$s[1],Number=1,Type=String,Description=\"Corresponding variant of the symbol\">";
my $a3 = "##FORMAT=<ID=SC,Number=1,Type=String,Description=\"Allele source\">";
my $time = scalar localtime();
my $a4 = "##source=$time $HEADER";
my $addinfo = "$a1\n$a2\n$a3\n$a4";

# read removed samples
my %remove = ();
if ($removesamp){
    my @aa = split /,/, $removesamp;
    foreach my $a (@aa){
        $remove{$a} = 1;
    }
}
if ($removefile){
    open RE, "<$removefile";
    while(<RE>){
        chomp;
        my @line = split /\t/, $_;
        $remove{$line[0]} = 1;
    }
}
################################################################ MAIN ###################################################################
#print STDOUT "$segregat\n";
my @rmpos = ();
my $prpos = ();
my $rmspos = ();
my %pos = ();
open IN, "<$input";
while(<IN>){
    chomp;
    my $in = $_;
    if (/^##/){
        print STDOUT "$in\n";
        next;
    }
    if (/^#/){
        print STDOUT "$addinfo\n";
        my @line = split /\t/, $in;
        my @samout = ();
        for (my $i = 9; $i <= $#line; $i ++){
            $pos{$line[$i]} = $i;
            if($remove{$line[$i]}){
                push @rmpos, $i;
                next;
            }
            next if ($line[$i] eq $p[0] or $line[$i] eq $p[1]);
            push @samout, $line[$i];
        }
        my $out1 = join "\t", @line[0..8];
        my $out2 = join "\t", @samout;
        print STDOUT "$out1\t$out2\n";
        push @rmpos, $pos{$p[0]};
        push @rmpos, $pos{$p[1]};
        $rmspos = join "\t", @rmpos;
        #rmpos cotains position of parents
        $prpos = "$pos{$p[0]}\t$pos{$p[1]}";
        #print STDOUT "$p[0]\t$pos{$p[0]}\t$rmspos\n";
        next;
    }
    my @allele_stat = split /\t/, &AlleleCount($rmspos, ,$prpos, $in);
    #print STDOUT "$in\n";
    #print STDOUT "$allele_stat[0]\n";
    next if($allele_stat[0] > $maxallele or $allele_stat[0] < $minallele);
    
    my @line = split /\t/, $in;
    #my @allele = split /,/, $line[4];
    #my $allele_num = $#allele + 2;
    #print STDOUT "$allele_num\n";
    #next if ($allele_num > $maxallele or $allele_num < $minallele);
    my $misspro = &MissCount ($rmspos, $in);
    next if ($misspro > $mmp);
    
    my $flag = "";
    if($method eq "m1"){
        $flag = &ParentsOnly($prpos, $in);
        if($flag eq "F"){
            next;
        }else{
            if ($maf){
                if ($allele_stat[1] >= $maf){
                    print STDOUT &ReGenotype($flag, $outsymbol, $rmspos, $in);
                }
            }else{
                print STDOUT &ReGenotype($flag, $outsymbol, $rmspos, $in);
            }
        }
    }elsif($method eq "m2"){
        $flag = &SingleParent($prpos, $in);
        #print STDOUT "\@$flag\n";
        next if ($flag eq "F");
        if ($maf){
            if ($allele_stat[1] >= $maf){
                #print STDOUT "$flag\n";
                print STDOUT &ReGenotype($flag, $outsymbol, $rmspos, $in);
            }
        }else{
            print STDOUT &ReGenotype($flag, $outsymbol, $rmspos, $in);
        }
    }
}

################################################################ SUB ###################################################################
sub MissCount {
    my ($rms, $in) = @_;
    my @line = split /\t/, $in;
    my @cc = split /\t/, $rms;
    my %rmsam = ();
    
    foreach my $c (@cc){
        $rmsam{$c} = 1;
    }
    my $miss = 0;
    my $sum = 0;
    for(my $i = 9; $i <= $#line; $i ++){
        next if ($rmsam{$i});
        if($line[$i] =~ /\.\/\./){
            $miss ++;
        }
        $sum ++;
    }
    my $out = $miss / $sum;
    return $out;
}
sub ParentsOnly {
    my ($ppos, $in) = @_;
    my @par = split /\t/, $ppos;
    my @line = split /\t/, $in;
    my %genotype = ();
    for (my $i = 9; $i <= $#line; $i ++){
        my @aa = split /:/, $line[$i];
        $genotype{$i} = $aa[0];
    }
    my $flag = 0;
    foreach my $p (@par){
        my @a = split /\//, $genotype{$p};
        if ($genotype{$p} =~ /\.\/\./){
            $flag = 1;
        }
        if ($a[0] ne $a[1]){
            $flag = 1;
        }
    }
    my $out = "";
    if ($flag == 1){
        $out = "F";
    }else{
        if ($genotype{$par[0]} eq $genotype{$par[1]}){
            $out = "F";
        }else{
            my $g1 = (split /\//, $genotype{$par[0]})[0];
            my $g2 = (split /\//, $genotype{$par[1]})[0];
            $out = "$g1\t$g2";
        }
    }
    return $out;
}

sub SingleParent {
    my ($ppos, $in) = @_;
    my @par = split /\t/, $ppos;
    my @line = split /\t/, $in;
    my %genotype = ();
    for (my $i = 9; $i <= $#line; $i ++){
        my @aa = split /:/, $line[$i];
        $genotype{$i} = $aa[0];
    }
    #remove het-site in parents
    my $flag = 0;
    foreach my $p (@par){
        my @a = split /\//, $genotype{$p};
        if ($a[0] ne $a[1]){
            $flag = 1;
        }
    }
    my $out = "F";
    if($flag == 0){
        if($genotype{$par[1]} eq "./." and $genotype{$par[0]} ne "./."){
            my $cc = "F";
            my $g1 = (split /\//, $genotype{$par[0]})[0];
            for (my $i = 9; $i <= $#line; $i ++){
                next if ($genotype{$i} =~ /\.\/\./);
                my @bb = split /\//, $genotype{$i};
                if ($bb[0] ne $g1){
                    $cc = $bb[0];
                }elsif($bb[1] ne $g1){
                    $cc = $bb[1];
                }
            }
            if ($cc eq "F"){
                $out = "F";
            }else{
                $out = "$g1\t$cc";
            }
        }elsif($genotype{$par[0]} eq "./." and $genotype{$par[1]} ne "./."){
            my $cc = "F";
            my $g2 = (split /\//, $genotype{$par[1]})[0];
            for (my $i = 9; $i <= $#line; $i ++){
                next if ($genotype{$i} eq "./.");
                my @bb = split /\//, $genotype{$i};
                if ($bb[0] ne $g2){
                    $cc = $bb[0];
                }elsif($bb[1] ne $g2){
                    $cc = $bb[1];
                }
            }
            if ($cc eq "F"){
                $out = "F";
            }else{
                $out = "$cc\t$g2";
            }
        }elsif($genotype{$par[0]} ne "./." and $genotype{$par[1]} ne "./."){
            my $flag = 0;
            foreach my $p (@par){
                my @gtpar = split /\//, $genotype{$p};
                if($gtpar[0] ne $gtpar[1]){
                    $flag = 1;
                }
            }
            if ($flag == 1){
                $out = "F";
            }else{
                if($genotype{$par[0]} ne $genotype{$par[1]}){
                    my $g1 = (split /\//, (split /:/, $line[$par[0]])[0])[0];
                    my $g2 = (split /\//, (split /:/, $line[$par[1]])[0])[0];
                    $out = "$g1\t$g2";
                }
            }
        }
    }
    return $out;
}


sub AlleleCount {
    my ($rms, $ppos, $in) = @_;
    my @par = split /\t/, $ppos;
    my @cc = split /\t/, $rms;
    my %rmsam = ();
    my %count = ();
    my %geno_count = ();
    my %gtype_count = ();
    my $sum = 0;
    my $het_sum = 0;
    foreach my $c (@cc){
        $rmsam{$c} = 1;
    }
    my @line = split /\t/, $in;
    my %genotype = ();
    for (my $i = 9; $i <= $#line; $i ++){
        my @aa = split /:/, $line[$i];
        $genotype{$i} = $aa[0];
    }
    my %gtall = ();
    foreach my $p (@par){
        my @a = split /\//, $genotype{$p};
        if ($genotype{$p} =~ /\.\/\./){
            next;
        }else{
            $gtall{$a[0]} ++;
            $gtall{$a[1]} ++;
        }
    }
   
    for (my $i = 9; $i <= $#line; $i ++){
        next if ($rmsam{$i});
        next if ($line[$i] =~ /\.\/\./);
        my @gtline = split /:/, $line[$i];
        my $ga = (split /\//, $gtline[0])[0];
        my $gb = (split /\//, $gtline[0])[1];
        $geno_count{$gtline[0]} ++;
        $gtype_count{$ga} ++;
        $gtype_count{$gb} ++;
        $gtall{$gb} ++;
        $gtall{$ga} ++;
        $sum ++;
        $het_sum ++ if ($ga ne $gb);
    }
    
    my $allele_num = keys %gtall;
    my $maf = 1;
    foreach my $a (sort keys %gtype_count){
        my $aa = $gtype_count{$a} / ($sum * 2);
        if ($aa <= $maf){
            $maf = $aa;
        }
    }
    my $out = "$allele_num\t$maf";
}


sub ReGenotype {
    my ($pgeno, $psym, $rms, $in) = @_;
    my @ps = split /,/, $psym;
    my @line = split /\t/, $in;
    my @cc = split /\t/, $rms;
    
    my %rmsam = ();
    foreach my $c (@cc){
        $rmsam{$c} = 1;
    }
    my $g1 = (split /\t/, $pgeno)[0];
    my $g2 = (split /\t/, $pgeno)[1];
    my %symbol = ();
    $symbol{$g1} = $ps[0];
    $symbol{$g2} = $ps[1];
    
    my %allele = ();
    $allele{0} = $line[3];
    my @ali = split /,/, $line[4];
    my $acount = 1;
    foreach my $ainfo (@ali){
        $allele{$acount} = $ainfo;
        $acount ++;
    }
    
    my @gtout = ();
    my $info = "$ps[0]=$allele{$g1};$ps[1]=$allele{$g2};";
    for (my $i = 9; $i <= $#line; $i ++){
        my $outgt = "";
        my @gtline = split /:/, $line[$i];
        if($rmsam{$i}){
            next;
        }else{
            if($gtline[0] eq "./."){
                $outgt = "./.";
            }else{
                my @gt = split /\//, $gtline[0];
                my @gto = ();
                push @gto, $symbol{$gt[0]};
                push @gto, $symbol{$gt[1]};
                @gto = sort @gto;
                $outgt = join "/", @gto;
            }
            push @gtout, $outgt;
        }
    }
    my $gt = join "\t", @gtout;
    my $site = join "\t", @line[0..6];
    my $out = "$site\t$info\tSC\t$gt\n";
    return $out;
}
