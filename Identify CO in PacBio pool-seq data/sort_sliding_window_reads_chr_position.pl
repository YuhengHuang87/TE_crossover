#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $file='';
my $j=0;my $ref="A6";my $que="A4";

my $outfile=$ref."_assembly_aln_benchmark_pool_posi_sorted";
open(OUT, ">$outfile");

my $outfile1=$ref."_assembly_aln_benchmark_pool_read_info_SNP_num";
open(OUT1, ">$outfile1");

my %A4_snp;my %A6_snp;my %A4_A6_convert;
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
    my $A4_posi = $d[1]."\t".$e[1];
    my $A6_posi = $b[0]."\t".$b[1];
		$A4_snp{$A4_posi}=$b[3]."\t".$b[4];
		$A6_snp{$A6_posi}=$b[3]."\t".$b[4];
		$A4_A6_convert{$A4_posi}=$A6_posi;
	}}}

my $read_ID; my $chr;my $read_name;
my @read_site; my @read_allele;
my %read; my $read_num=0; my $dif_num=0; my $rep_num=0; my %dif_al; my $check=0;my $read_info_snp=0; my $dif_chr=0;my $same_chr=0;
$file = "aln_".$ref."_benchmark_pool_sorted_que20_read_posi";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		my $site=$b[1]."\t".$b[4];
		my $snp_position='';
		if (exists $A6_snp{$site}){
			$snp_position=$b[0]."\t".$site;
		}elsif(exists $A4_snp{$site}){
			$snp_position=$b[0]."\t".$A4_A6_convert{$site};
		}
		if ($snp_position ne ''){
			my $allele=$snp_position."\t".$b[1]."\t".$b[2]."\t".$b[3]."\t".$b[4]."\t".$b[5]."\t".$b[6]."\t".$b[7];
		if ($b[0] eq $read_ID){
			if (exists $read{$snp_position}){
				$rep_num++;
				if ($b[5] ne $read{$snp_position}){
					#print "$count\t","$read{$read_position}\n";
					$dif_num++;
					$dif_al{$snp_position}=$b[3];
				}
			}else{
				my @loca=split("\t", $snp_position);
				if ($loca[1] eq $chr){
					$same_chr++;
				$read_info_snp++;
				$read{$snp_position}=$b[5];
				push @read_site, $loca[2];
				push @read_allele, $allele;
			}else{
				$dif_chr++;
			}
		}
		}else{
			print OUT1 "$read_name\t","$read_info_snp\n";
			$read_info_snp=1;
			my @idx = sort { $read_site[$a] <=> $read_site[$b] } 0 .. $#read_site;

			@read_site = @read_site[@idx];
			@read_allele = @read_allele[@idx];
			for (my $i=0; $i < @read_site; $i++) {
				my $site=$read_name."\t".$read_site[$i];
				if (exists $dif_al{$site}){
					$check++;
					#print OUT "$read_allele[$i]\t","dup\n";
				}else{
				print OUT "$read_allele[$i]\t","unique\n";
			}
		}
			my @f=split("\t", $snp_position);
			$read_ID=$f[0];$chr=$f[1];$read_name=$read_ID."\t".$chr;
			@read_site=($b[4]);@read_allele=($allele);$read_num++;
		}
	}}
print "$read_num\t","$rep_num\t","$dif_num\t","$check\t","$same_chr\t","$dif_chr\n";

my @snp_leng=("4","10");
foreach my $snp_num (@snp_leng){
my $outfile3=$ref."_assembly_region_sliding_window_benchmark_pool_prop_posi_".$snp_num.".txt";
open(OUT3, ">$outfile3");

my @win_allele=();my @array_posi=();
my $read_name=''; my $site_num=0; my $real_allele; my $win_perce; my $heter_perce; my $posi;my $region; my $last_region;
my $win_snp=0; my $file='';

$file = $ref."_assembly_aln_benchmark_pool_posi_sorted";
#$ref."_assembly_aln_benchmark_pool_posi_sorted";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
  $posi = $b[0]."\t".$b[1]."\t".$b[2]."\t".$b[4]."\t".$b[5];
	$real_allele=$b[7];
  $region=$b[0]."\t".$b[1]."\t".$b[4];

	#if ($b[0] eq $read_name){
    if ($region eq $last_region){
#my @allele=split("\t", $snp{$posi});
	my @allele=($b[8],$b[9]);
	if ($real_allele eq $allele[0]){
  push @win_allele,1;
	push @array_posi,$b[2];
  $site_num++;
}elsif ($real_allele eq $allele[1]){
	push @win_allele,0;
	push @array_posi,$b[2];
  $site_num++;
}

if (@win_allele  == $snp_num){
$win_perce=sum(@win_allele)/$snp_num;
print OUT "$posi\t","$site_num\t","$win_perce\t","$array_posi[0]\t","$array_posi[-1]\n";
splice @win_allele,0,1;
splice @array_posi,0,1;
}
}else{
	@win_allele=();@array_posi=();
	my @allele=($b[8],$b[9]);
		if ($real_allele eq $allele[0]){
	  push @win_allele,1;
		push @array_posi,$b[2];
    $site_num=1;
	}elsif ($real_allele eq $allele[1]){
		push @win_allele,0;
		push @array_posi,$b[2];
    $site_num=1;
	}
	#$read_name=$b[0];
  $last_region=$region;
}}
close I;
}