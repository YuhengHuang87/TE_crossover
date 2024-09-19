#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $file = '';

my $chromo="CM010569.1";
my $alter_chr="CM010541.1";
my $minimum= 144265;
my $maximum = 18689904;
#my @euch=("CM010569.1_144265_18689904","CM010570.1_251300_19760175","CM010571.1_8030989_23829870","CM010572.1_164091_18449201","CM010573.1_8729782_31194759");

my $indiv=400;
my $loca;

my $outfile1="A6_A4_recombine_arms_each_".$chromo.".fasta";
open(OUT1, ">$outfile1");

my $outfile="A6_A4_recombine_individual_arms_".$chromo."_location";
open(OUT, ">$outfile");

my $nuc_A6;
my $j;
$file = "/dfs7/grylee/yuhenh3/recombination_rate/A6_genomic.fna";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
if($count =~ m/>/){
my @a = split(" ", $count);
my @b = split(">", $a[0]);
print "$b[1]\n";
if($b[1] eq $chromo){
$nuc_A6="";
$j=1;
}else{
$j=0;
}
}else{
  if (($j==1)){
$nuc_A6=$nuc_A6.$count;
}}}

my $nuc_A4;
$file = "/dfs7/grylee/yuhenh3/recombination_rate/A4_genomic.fna";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
		if($count =~ m/>/){
		my @a = split(" ", $count);
		my @b = split(">", $a[0]);
		if($b[1] eq $alter_chr){
		$nuc_A4="";
		$j=1;
		}else{
		$j=0;
		}
		}else{
		  if (($j==1)){
		$nuc_A4=$nuc_A4.$count;
		}}}

my @chr_nuc_A6=split("", $nuc_A6);
my @chr_nuc_A4=split("", $nuc_A4);

my $i=1;
while ($i<=$indiv){
$loca = $minimum + int(rand($maximum - $minimum));

my $ref_str=''; my $ref_end='';
my $k=0; my $alter_chro='';my $alter_loca='';
my @allele_ref; my @allele_alt;
$file="/dfs7/grylee/yuhenh3/recombination_rate/A6_refer_A4_query_asm10_update.var.vcf";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
      if ($b[0] eq $chromo){
				if ($k==0){
					$ref_str=$b[1];
				}else{
					$ref_end=$b[1];

					if (($loca > $ref_str)&&($loca < $ref_end)){
						my $dis1= scalar (@allele_ref); my $dis2 = scalar (@allele_alt);
						if (($loca > $ref_str+$dis1-1)&&($loca > $ref_str+$dis2-1)){
						my @a = split(";", $b[-3]);
						my @c = split("=", $a[0]);
						my @d = split("=", $a[1]);
						$alter_chro = $c[1];
						$alter_loca = $d[1]+$loca-$ref_str;
					}else{
						print "$count\t","$loca\n";
					}
        }}
		@allele_ref = split("", $b[3]);@allele_alt = split("", $b[4]);
				$ref_str=$b[1];
				$k=1;
			}
		}

if ($alter_loca ne ''){
print OUT1 ">","$chromo","_","$i","_","$loca\n";
print OUT "$i\t","$chromo\t","$loca\t","$alter_chro\t","$alter_loca\n";

my $is_odd = $loca % 2;
if ($is_odd==0){
for (my $m=0; $m<$loca; $m+=1){
print OUT1 "$chr_nuc_A6[$m]";
}
for (my $n=$alter_loca; $n<@chr_nuc_A4; $n+=1){
print OUT1 "$chr_nuc_A4[$n]";
}
print OUT1 "\n";
}else{
	for (my $m=0; $m<$alter_loca; $m+=1){
	print OUT1 "$chr_nuc_A4[$m]";
	}
	for (my $n=$loca; $n<@chr_nuc_A6; $n+=1){
	print OUT1 "$chr_nuc_A6[$n]";
	}
print OUT1 "\n";
}
$i++;
}}
