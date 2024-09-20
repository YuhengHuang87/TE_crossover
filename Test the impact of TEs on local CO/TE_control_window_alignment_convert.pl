#!/usr/bin/perl
use strict;
use warnings;

my @TE_ref=();
my $a=0; my $t=0;
my $homolog_location="A6";
my $strain_focal="A7";
my $win_size=5000;my $control_win_size=$win_size*2;
#my $outfile=$strain_focal."focal_5000bp_distance_".$homolog_location."_location_TE_".$win_size.".txt";
my $outfile=$strain_focal."_focal_control_alter_".$homolog_location."_location_window_".$control_win_size.".txt";
open(OUT, ">$outfile");


#my $infile1 = $strain_focal."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic.txt";
my $infile1=$strain_focal."_crossover_num_control_included_window_exclude_TE_overlapped_average_".$control_win_size.".txt";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
#my $locat=$b[0]."\t".$b[1]."\t".$b[2]."\t".$b[3]; # for TE windows
my $locat=$b[0]."\t".$b[0]."\t".$b[1]."\t".$b[2]; #for control windows
push @TE_ref, $locat;
$t++;
}

for (my $k=0; $k<@TE_ref; $k++){
my @posi = split("\t", $TE_ref[$k]);
my $chr=$posi[1]; my $left=$posi[2];my $right=$posi[3];
my $i=0; my $j=0;
my $ref_str=''; my $ref_end=''; my $que_str=''; my $que_end='';
my $te_q_s; my $te_q_e;

my $infile2=$strain_focal."_refer_".$homolog_location."_query_asm10_update.var.txt"; 
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
if($b[0] eq "V"){
if ($b[1] eq $chr){
$ref_end=$b[2]; #need to change
$que_end=$b[9];

#$ref_end=$b[9];
#$que_end=$b[2];

if ($ref_str ne ""){
if(($left>$ref_str)&&($left<$ref_end)){
my $up_dis_s=$left-$ref_str;
$te_q_s=$que_str+$up_dis_s;
$i++;
}elsif(($left>=$b[2])&&($left<=$b[3])){
	$te_q_s=$b[-3];
	$i++;
}
if(($right>$ref_str)&&($right<$ref_end)){
my $up_dis_e=$right-$ref_str;
$te_q_e=$que_str+$up_dis_e;
$j++;
}elsif(($right>=$b[2])&&($right<=$b[3])){
	$te_q_e=$b[-2];
	$j++;
}
#if (($i==1)&&($j==1)){
	#if ($i==1){
	if (($i>0)&&($j>0)){
print OUT "$TE_ref[$k]\t";
print OUT "$b[8]\t","$te_q_s\t","$te_q_e\n";
#print OUT "$te_q_s\n";
$a++;
last;
}}
$ref_str=$b[3];#need to change
$que_str=$b[10];

#$ref_str=$b[10];
#$que_str=$b[3];

}}}}
print "$a\t","$t\n";
