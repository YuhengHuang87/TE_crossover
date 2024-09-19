#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
my $hmd_unit=25; my $win_leng=1000; my $num_measure = 10; my $win_con=2; my $fold_dif=1;
my @name=("A6","A7");
my @combine=('1_2_4','2_3_4_5');
my $strain="A6"; my $alter="A7";

my $index; my $rep;
$index=1;#need to match strain or alter
$rep=0;

my @index_pair=split("_", $combine[$index]);
my $p= $index_pair[$rep];

my $i=0;my $j=0;
#my $outfile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/".$name[$index]."_".$p."_assigned_TE_spreading_HMD_25_".$win_leng."_".$num_measure."_measures_median_include_ambiguous.txt";
my $outfile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/focall_".$strain."_alter_".$alter."_".$p."_TE_spreading_HMD_25_".$win_leng."_".$num_measure."_include_ambiguous.txt";
open(OUT, ">$outfile1");

my %background;my %hmd;
#my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_mean_control".$name[$index]."_".$p."_exclu_within_5000bp_distance_nearby_TE_euchromatic.txt";
my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_mean_focal_".$strain."_alter_".$alter."_".$p.".txt";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $te_posi=$b[0]."\t".$b[1]."\t".$b[2];
if (($b[3] > 0)&&($b[-1] >= 1000)){
$background{$te_posi}=$b[3];
}}

my $infile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/".$name[$index]."_".$p."_25bp_avgCov_unique_mapped";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $locat=$b[0]."\t".$b[1];
$i++;
$hmd{$locat}=$b[2];
}

my $con_mean;my $k=0;

#my $infile2 = "/dfs7/grylee/yuhenh3/recombination_rate/".$name[$index]."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic_blast_euchromatic_include_ambiguous_10";
my $infile2="/dfs7/grylee/yuhenh3/recombination_rate/".$strain."focal_5000bp_distance_".$alter."_location_TE_5000.txt";

open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		
my $chr=$b[-3]; my $star=$b[-2]; my $end=$b[-1];
my $te_posi=$chr."\t".$star."\t".$end;

if (exists $background{$te_posi}){
$con_mean=$background{$te_posi};
	my %dis_l_enrich; my %dis_r_enrich;my $win; my $dis; my $dis_raw;my @TE_wi_hmd=();
	foreach my $key (keys %hmd) {
my @posi = split("\t", $key);
if ($posi[0] eq $chr){
$win=$posi[1];
if (($win+$hmd_unit >= $star - 50000)&&($win <= $end+50000)){
if ($win+$hmd_unit < $star){
$dis = int(($win+$hmd_unit-$star)/$win_leng);
if(exists $dis_l_enrich{$dis}){
$dis_l_enrich{$dis}=$dis_l_enrich{$dis}."\t".$hmd{$key};
}else{
$dis_l_enrich{$dis}=$hmd{$key};
}
}elsif ($win > $end){
$dis = int(($win-$end)/$win_leng);
if(exists $dis_r_enrich{$dis}){
$dis_r_enrich{$dis}=$dis_r_enrich{$dis}."\t".$hmd{$key};
}else{
$dis_r_enrich{$dis}=$hmd{$key};
}
}else{
push @TE_wi_hmd, $hmd{$key};
}
}}}


my $left_ext=''; my $l_count=0; my $l_miss=0;
foreach my $key (sort {$b <=> $a } keys %dis_l_enrich){
my @near_w = split("\t", $dis_l_enrich{$key});
if (@near_w>=$num_measure){
if (mean(@near_w)<=$con_mean*$fold_dif){
$l_count++;
}else{
$l_count=0;
}
}else{
$l_miss++;
}

if ($l_count==$win_con){
$left_ext=($key-$l_miss)*$win_leng;
last;
}}

my $right_ext='';my $r_count=0;my $r_miss=0;
foreach my $key (sort {$a <=> $b } keys %dis_r_enrich){
my @near_w = split("\t", $dis_r_enrich{$key});
if (@near_w>=$num_measure){
if (mean(@near_w)<=$con_mean*$fold_dif){
$r_count++;
}else{
$r_count=0;
}
}else{
	$r_miss++;
}

if ($r_count==$win_con){
$right_ext=($key-$r_miss)*$win_leng;
last;
}}

my $final_ext;my $na=0;
if (($left_ext ne '')&&($right_ext ne '')){
$final_ext = (abs($left_ext)+abs($right_ext))/2;
}elsif(($left_ext eq '')&&($right_ext ne '')){
$final_ext=abs($right_ext);$left_ext='NA';
}elsif(($left_ext ne '')&&($right_ext eq '')){
$final_ext=abs($left_ext);$right_ext='NA';
}else{
$na++;
}
if ($na==0){

if (@TE_wi_hmd>0){
my $avg_TE_hmd=mean(@TE_wi_hmd);
my $median_TE_hmd=median(@TE_wi_hmd);
print OUT "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\t","$left_ext\t","$right_ext\t","$final_ext\t","$l_miss\t","$r_miss\t","$avg_TE_hmd\t","$median_TE_hmd\n";
}else{
print OUT "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\t","$left_ext\t","$right_ext\t","$final_ext\t","$l_miss\t","$r_miss\t","NA\t","NA\n";
}}
}
}
#}
#}

my $outfile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/Mass_".$strain."_alter_".$alter."_".$p."_assigned_TE_height_HMD_25_".$win_leng."_".$num_measure."_measures_mean_include_ambiguous.txt";
open(OUT1, ">$outfile1");

my %background;my %hmd;
#my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_mean_control".$strain."_".$p."_exclu_within_5000bp_distance_nearby_TE_euchromatic.txt";
my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_mean_focal_".$strain."_alter_".$alter."_".$p.".txt";

#####
#my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_mean_focal_".$name."_alter_".$alter."_".$p."_crossover_num_TE_included.txt";

#my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_control_focal_".$name."_".$p."_".$control_win_size.".txt";
#my $infile="/dfs7/grylee/yuhenh3/CUT_tag_experiment/HMD_unique_mapped_local_control_alter_".$alter."_".$p."_".$control_win_size.".txt";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $te_posi=$b[0]."\t".$b[1]."\t".$b[2];
if (($b[3] > 0)&&($b[-1] >= 1000)){
$background{$te_posi}=$b[3];
}}

my $infile1="/dfs7/grylee/yuhenh3/CUT_tag_experiment/".$name[$index]."_".$p."_25bp_avgCov_unique_mapped";#always this file but change the name
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $locat=$b[0]."\t".$b[1];
$hmd{$locat}=$b[2];
}

my $median;my $k=0;
#my $infile2 ="/dfs7/grylee/yuhenh3/CUT_tag_experiment/".$strain."_".$p."_assigned_TE_spreading_HMD_25_".$win_leng."_".$num_measure."_measures_median_include_ambiguous.txt";
my $infile2="/dfs7/grylee/yuhenh3/CUT_tag_experiment/focall_".$strain."_alter_".$alter."_".$p."_TE_spreading_HMD_25_".$win_leng."_".$num_measure."_include_ambiguous.txt";

#########
#my $infile2 = "/dfs7/grylee/yuhenh3/recombination_rate/".$name."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic_blast_euchromatic_include_ambiguous_10";
########
#my $infile2="/dfs7/grylee/yuhenh3/recombination_rate/focal_".$name."_alter_".$alter."_crossover_num_TE_included_window_exclude_TE_overlapped_average_".$win_size.".txt";
########
#my $infile2="/dfs7/grylee/yuhenh3/recombination_rate/focal_".$name."_alter_".$alter."_crossover_num_control_included_window_exclude_TE_overlapped_average_".$control_win_size.".txt";
#######

open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
#my $chr=$b[1]; my $star=$b[2]; my $end=$b[3]; #for focal TEs
my $chr=$b[4]; my $star=$b[5]; my $end=$b[6]; #for alternative homologs

my $left_ext = $b[7]; my $right_ext = $b[8];
my $left_step = $left_ext/1000; my $right_step = $right_ext/1000;

my $te_posi=$chr."\t".$star."\t".$end;
	my $l_mass=0; my $r_mass=0;
    my $win; 
	if (exists $background{$te_posi}){
my $mean_control=$background{$te_posi};
    for (my $i = 0; $i < abs($left_step); $i++){
        my @K9_enrich=();
	foreach my $key (keys %hmd) {
    my @posi = split("\t", $key);
    if ($posi[0] eq $chr){
        $win=$posi[1];
        if (($win >= $star + $left_ext + $i*1000)&&($win+$hmd_unit < $star + $left_ext + $i*1000 + 1000)){
            push @K9_enrich, $hmd{$key}/$mean_control;
        }
    }}
    if (@K9_enrich>0){
    my $K9_mean = mean(@K9_enrich); my $K9_median = median(@K9_enrich);
    $l_mass=$l_mass+$K9_mean;
    }}
    

for (my $i = 0; $i < abs($right_step); $i++){
        my @K9_enrich=();
	foreach my $key (keys %hmd) {
    my @posi = split("\t", $key);
    if ($posi[0] eq $chr){
        $win=$posi[1];
        if (($win >= $end + $i*1000)&&($win+$hmd_unit < $end + $i*1000 + 1000)){
            push @K9_enrich, $hmd{$key}/$mean_control;
        }
    }}
        if (@K9_enrich>0){
    my $K9_mean = mean(@K9_enrich); my $K9_median = median(@K9_enrich);
    $r_mass=$r_mass+$K9_mean;
    }}
print OUT1 "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\t","$b[7]\t","$b[8]\t","$l_mass\t","$r_mass\n";
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
