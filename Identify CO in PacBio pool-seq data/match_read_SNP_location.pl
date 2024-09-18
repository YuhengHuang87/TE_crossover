#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

#please check the format of CIGAR in sam files and adjust the code from line 78-91

my $k;
my $file = '';
$k="00";

my $num=0; my $alter=0;
my $line=0; my %name;
my $ref="A6";my $que="A4";
my $outfile=$k."_aln_".$ref."_benchmark_pool_sorted_que20_read_posi";
open(OUT, ">$outfile");

my $num_site=0;my %A6_snp;my %A4_snp;
$file = $ref."_refer_".$que."_query_asm10_update_SNP_shared_assem_Illumia_two.vcf";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    if($count =~ m/#/){
}else{
    my @b = split("\t", $count);
    my @ref=split("", $b[3]);    my @alt=split("", $b[4]);
    if ((@ref == 1)&&(@alt == 1)){
      my @c=split(";", $b[7]);
      my @d= split("=", $c[0]);
      my @e= split("=", $c[1]);
#print "$c[1]\n";
    my $A4_posi = $d[1]."\t".$e[1];
    my $A6_posi = $b[0]."\t".$b[1];
    my $nuclotide=$b[3]."\t".$b[4];
    $A6_snp{$A6_posi}=$nuclotide;
    #$A4_snp{$A4_posi}=$nuclotide;
    $num_site++;
}}}
#print "$num_site\n";

my %read;my %num_dup; my $dup=0; my $total_mapped=0; my $SA=0;

$file ="aln_A7_A7Pool_batch3_sorted_que20_prefix".$k;
open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
my @a = split("\t", $count);
if ($a[2] =~ m/CM0105/){
my $start=$a[3];
my @b=split("",$a[5]);
my @read_list;my $len=0;
my @rep; my @type;
my $num=""; my $index=0;
foreach my $b (@b){
	if ($b =~ m/[0-9]/){
$num=$num.$b;
}else{
$rep[$index]=$num;
$type[$index]=$b;
$index++;
$num="";
}}


for (my $i=0; $i < scalar(@rep); $i++) {
@read_list=(@read_list,($type[$i])x$rep[$i]);
if ($type[$i] =~ m/=|X|D|N/){
  #if ($type[$i] =~ m/M|D|N/){
$len=$len+$rep[$i];
}}

my @nucl_list=split("",$a[9]);
my $nucl_index=-1;
my $nucl_position=$start-1;
my $read_position=0;

foreach my $base_type (@read_list){
if($base_type =~ m/=|X|D|N/){
#if($base_type =~ m/M|D|N/){
$nucl_position++;
}
if($base_type =~ m/=|X|I|S/){
#if($base_type =~ m/M|I|S/){
$nucl_index++;
}
if($base_type =~ m/=|X|I|S|H/){
#if($base_type =~ m/M|I|S|H/){
$read_position++;
}
if(($base_type eq "=")||($base_type eq "X")){
#if(($base_type eq "M")){
my $read_allele=$a[0]."\t".$a[2]."\t".$a[3]."\t".$read_position."\t".$nucl_position."\t".$nucl_list[$nucl_index];
my $site=$a[2]."\t".$nucl_position;
	if($a[2] =~ m/105/){
if (exists $A6_snp{$site}){
print OUT "$read_allele\t","$A6_snp{$site}\n";
}

}
}
}
}
}
