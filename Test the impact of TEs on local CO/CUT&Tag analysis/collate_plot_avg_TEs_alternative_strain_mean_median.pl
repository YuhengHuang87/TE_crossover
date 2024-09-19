#!/usr/bin/perl
use strict;
use warnings;

my @name=("A6","A7");
my @combine=('1_2_4','2_3_4');

my $index=0; my $alter=1;
my $rep=0;
my $step=1000; my $TE_exclude_step = 10000;
#my @index_pair=split("_", $combine[$index]);
my @index_pair=split("_", $combine[$alter]);

my $p= $index_pair[$rep];

my $i=0;my $j=0;
#my $outfile1="K9_Collate_betw_strain_".$name[$index]."_".$p."_TE_200bp_50kb_win_".$step."_".$TE_exclude_step."_local_K9_standardized_mean_include_ambiguous.txt";#change the name
my $outfile1="K9_Collate_betw_strain_".$name[$index]."_alternative_".$name[$alter]."_".$p."_TE_200bp_50kb_win_".$step."_".$TE_exclude_step."_local_K9_standardized_mean_include_ambiguous.txt";

open(OUT, ">$outfile1");

my %background;my %hmd;
#my $infile="HMD_unique_mapped_local_mean_control".$name[$index]."_".$p."_exclu_within_5000bp_distance_nearby_TE_euchromatic.txt";#change the name
my $infile="HMD_unique_mapped_local_mean_control_".$name[$index]."_focal_".$name[$alter]."_".$p."_exclu_within_5000bp_distance_nearby_TE_euchromatic.txt";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $te_posi=$b[0]."\t".$b[1]."\t".$b[2];
if (($b[3] > 0)&&($b[-1] >= 1000)){
$background{$te_posi}=$b[3];
}}

#my $infile1=$name[$index]."_".$p."_25bp_avgCov_unique_mapped";
my $infile1=$name[$alter]."_".$p."_25bp_avgCov_unique_mapped";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $locat=$b[0]."\t".$b[1];
$i++;
$hmd{$locat}=$b[2];
}


my @alter_TE; my $loca;
my $file =$name[$alter]."_focal_".$name[$index]."_location_TE_assigned_include_ambiguous.txt";
open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    $loca = $a[1]."\t".$a[2]."\t".$a[3];
    push @alter_TE, $loca;
}

my @focal_loca; 
my @alter_loca; my $alter_loca; my $overlap=0;
$file =$name[$index]."_focal_".$name[$alter]."_location_TE_assigned_include_ambiguous.txt";
open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    $overlap=0;
    my @a = split("\t", $count);
    $loca = $a[1]."\t".$a[2]."\t".$a[3];
    $alter_loca = $a[-3]."\t".$a[-2]."\t".$a[-1];
    my $leng = $a[3]-$a[2]+1;
    if ($leng>=200){
    for (my $j=0; $j<@alter_TE; $j++){
my @posi = split("\t", $alter_TE[$j]);
my $TE_left=$posi[1]; my $TE_right=$posi[2];
if ($posi[0] eq $a[-3]){
if (($TE_left > $a[-1]+$TE_exclude_step)||($TE_right < $a[-2]-$TE_exclude_step)){
}else{
    $overlap++;
}
}}
   if ($overlap==0){
        push @focal_loca, $loca;
        push @alter_loca, $alter_loca;
  }
}
}


my %sum_enrich; my %win_num; my $win; my $dis; my $dis_raw;my $median;
    for (my $k=0; $k<@alter_loca; $k++){
		my @b = split("\t", $alter_loca[$k]);

my $te_posi=$b[0]."\t".$b[1]."\t".$b[2];
my $left=$b[1]; my $right=$b[2];my $chr=$b[0];
#print "$alter_loca[$k]\n";
if (exists $background{$te_posi}){
my %K9_TE; 
	foreach my $key (keys %hmd) {
my $normalized_K9=$hmd{$key}/$background{$te_posi};
my @posi = split("\t", $key);
if ($posi[0] eq $chr){
$win=$posi[1]; my $inter=($left+$right)/2;
if (($win >= $left - 40000)&&($win <= $right+40000)){
if ($win < $left){
$dis_raw = $win-$left;
}
if ($win > $right){
$dis_raw = $win-$right;
}
$dis=int($dis_raw/$step);
if(exists $K9_TE{$dis}){
$K9_TE{$dis}=$K9_TE{$dis}."\t".$normalized_K9;
}else{
$K9_TE{$dis}=$normalized_K9;
}}}
}

foreach my $key (sort {$a <=> $b } keys %K9_TE){
my @enrich=split("\t", $K9_TE{$key});
my $mean_enrich=mean(@enrich);
if(exists $sum_enrich{$key}){
$sum_enrich{$key}=$sum_enrich{$key}."\t".$mean_enrich;
$win_num{$key}++;
}else{
$sum_enrich{$key}=$mean_enrich;
$win_num{$key}=1;
}
}
}}

foreach my $key (sort {$a <=> $b } keys %sum_enrich){
my @enrich=split("\t", $sum_enrich{$key});
my $med_enrich=median(@enrich);
my $mean_enrich=mean(@enrich);
if ($win_num{$key}>200){
print OUT "$key\t","$mean_enrich\t","$med_enrich\t","$win_num{$key}\n";
}
}



sub mean {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
