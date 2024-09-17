#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my $snp_num=50;
my $edge_bp=200000; my $edge_snp=50;

my $file=''; my %snp;
my @ID=(2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24,27,28,29,30,31,33,34,35,37,38,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,63,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,90,91,92,93,94,95,96,97,98,99,100,101,103,104,105,106,107,108,110,111,113,114,115,116,117,118,119,120,121,122,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,164,165,166,168,169,170,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,198,199,200,201,202,203,204,205,206,207,209,210,211,212,213,214,215);

$file = "/dfs3b/grylee/yuhenh3/recombination_rate/A6_refer_A4_query_asm10_update_SNP_shared_assem_Illumia_two.vcf";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    if($count =~ m/#/){
}else{
    my @b = split("\t", $count);
    my $site = $b[0]."\t".$b[1];
		my $allele=$b[3]."\t".$b[4];
    $snp{$site}=$allele;
}
}

my $outfile="A4_A6_assembly_Illumina_geno_sliding_window_A6_breakpoint_all_Illumia_update_".$snp_num."_".$edge_bp."_".$edge_snp;
open(OUT, ">$outfile");

foreach my $ID (@ID){
  my $outfile1=$ID."_A6_parenatal_A4_A6_Homo_Heter_all_sliding_win_".$snp_num."SNPs_percent_shared_Illumia_update.txt";
	open(OUT, ">$outfile1");

my @win_allele=();my @array_posi=();
my $chr=''; my $crossover_read=0; my $total_num=0; my $win_perce; my $heter_perce; my $posi;
my $win_snp=0;

$file = $ID."_A6_genomic_sorted.g.vcf";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
	chomp($count);
	if($count =~ m/#/){
}else{
	my @b = split("\t", $count);
  $posi = $b[0]."\t".$b[1];
		my @ref=split("", $b[3]);    my @alt=split("", $b[4]);
		if ((@ref == 1)&&(@alt == 1)){
			my @sam=split(":", $b[-1]);
			#if ($sam[2] >=10){
			my @info=split(";", $b[-3]);
			my @geno=split("=", $info[1]);
			if (exists $snp{$posi}){
				my @al=split("\t", $snp{$posi});
				my $A6_allele=$al[0];my $A4_allele=$al[1];
				if (($b[3] eq $A6_allele)&&($b[4] eq $A4_allele)){
			if ($b[0] eq $chr){
			if ($geno[1] == 0.5){
				$total_num++;
				push @win_allele,1;
				push @array_posi,$b[1];
			}elsif($geno[1] == 1){
				$total_num++;
				push @win_allele,0;
				push @array_posi,$b[1];
				}

if (@win_allele  == $snp_num){
$win_perce=sum(@win_allele)/$snp_num;
print OUT "$posi\t","$win_perce\t","$array_posi[0]\t","$array_posi[-1]\n";
splice @win_allele,0,1;
splice @array_posi,0,1;
}
}else{
	@win_allele=();@array_posi=();
	if ($geno[1] == 0.5){
		push @win_allele,1;
		push @array_posi,$b[1];
	}elsif($geno[1] == 1){
		push @win_allele,0;
		push @array_posi,$b[1];
		}
	$chr=$b[0];
	}
}}
#}
}
}}
close I;
print "$ID\t","$total_num\n";


my $num_read=0; my $dup_read=0;
my %candid_read; my @loca=(); my $posi;

$file =$ID."_A6_parenatal_A4_A6_Homo_Heter_all_sliding_win_".$snp_num."SNPs_percent_shared_Illumia_update.txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		$posi = $b[0]."\t".$b[1]."\t".$b[3]."\t".$b[4];
    if ($b[2] == 0.5){
			if (exists $candid_read{$b[0]}){
				push @loca, $posi;
				$dup_read++;
			}else{
			push @loca, $posi;
		$candid_read{$b[0]}=1;
		$num_read++;
}}}
print "$ID\t","$dup_read\t","$num_read\n";

my $left_count=0; my $right_count=0; my $left_consist=0;my $right_consist=0;
for (my $i=0; $i<@loca; $i+=1){
	my @a = split("\t", $loca[$i]);
	$left_count=0; $right_count=0; $left_consist=0;$right_consist=0;

$file = $ID."_A6_parenatal_A4_A6_Homo_Heter_all_sliding_win_".$snp_num."SNPs_percent_shared_Illumia_update.txt";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
	if (exists $candid_read{$b[0]}){
		if ($b[0] eq $a[0]){
		if (($b[4] < $a[3])&&($b[4] > $a[3]-$edge_bp)){
			if ($b[2] < 0.5){
				$left_consist=$left_consist-1;
				$left_count++;
			}elsif($b[2] > 0.5){
				$left_consist=$left_consist+1;
				$left_count++;
			}
}elsif (($b[4] > $a[3])&&($b[4] < $a[3]+$edge_bp)){
	if ($b[2] < 0.5){
		$right_consist=$right_consist-1;
		$right_count++;
	}elsif($b[2] > 0.5){
		$right_consist=$right_consist+1;
		$right_count++;
	}
}}}}
close I;
#print "$loca[$i]\t","$left_count\t","$right_count\t","$left_consist\t","$right_consist\n";
if (($left_count >= $edge_snp)&&($right_count >= $edge_snp)){
	if ((abs($left_consist) == $left_count)&&(abs($right_consist) == $right_count)){
		if ($left_consist*$right_consist<0){
		print OUT "$ID\t","$loca[$i]\t","$left_consist\t","$right_consist\n";
	}}
}
}
}
