#!/usr/bin/perl
use strict;
use warnings;
my $snp_num=10;
my $distance = 5000;
my $ref="A6";my $que="A4"; 
my $pool="benchmark"; #change it to A6, A7 experiment pool
my $read_name=''; my $pass_num=0; my $read_num=0;
my $win_snp=0; my $start=0;my $end=0;
my $file='';my %ID; 

$file = $ref."_assembly_aln_".$pool."_pool_posi_sorted";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
if ($b[0] eq $read_name){
    $win_snp++;
    $end = $b[2];
}else{
    if ($win_snp >$snp_num){
        if (($end-$start)>$distance){
        $ID{$read_name}=1;
        $pass_num++;
    }}
    $read_num++;
    $read_name = $b[0];
    $win_snp=0;
    $start = $b[2];
    $end = $b[2];
}
}
close I;


my $outfile=$ref."_assembly_".$pool."_pool_read_pass_snp_".$snp_num"._depth.txt";
open(OUT, ">$outfile");

my $i=0;
my %snp_cov;
my $infile2=$ref."_assembly_aln_".$pool."_pool_posi_sorted";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
my @b = split("\t", $count);
my $id = $b[0];
#my $id = $b[0]."\t".$b[1];
if (exists $ID{$id}){
my $posi = $b[1]."\t".$b[2];
if (exists $snp_cov{$posi}){
    $snp_cov{$posi}=$snp_cov{$posi}+1;
}else{
    $snp_cov{$posi}=1;
}
}
}

my $file = $ref."_refer_".$que."_query_asm10_update_SNP_shared_assem_Illumia_two.vcf";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
    my $A6_posi = $b[0]."\t".$b[1];
    my $cov;
    if (exists $snp_cov{$A6_posi}){
        $cov = $snp_cov{$A6_posi};
    }else{
        $cov = 0;
    }
print OUT "$A6_posi\t","$cov\n";
}