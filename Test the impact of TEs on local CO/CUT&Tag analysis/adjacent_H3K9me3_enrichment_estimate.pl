#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
my $hmd_unit=25; my $win_leng=5000; my $num_measure=50;
my $win_size=5000;  my $dis1=30000;my $dis2=40000;
my @name=("A6","A7");
my @combine=('1_2_4','2_3_4');
my $strain="A6"; my $alter="A7";#need to match the index below

my $index; my $rep;
$index=1;#need to match focal strain or alter depending on what it's running
$rep=0;

my @index_pair=split("_", $combine[$index]);
my $p= $index_pair[$rep];
my $outfile1="nearest_distance_focal_".$strain."_".$p."_alter_".$alter."_height_HMD_25_".$win_leng."_".$num_measure."_measures_mean.txt"; 
open(OUT, ">$outfile1");

my %background;my %hmd;

my $infile="HMD_unique_mapped_local_nearest_distance_alter_".$alter."_".$p.".txt";
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $te_posi=$b[0]."\t".$b[1]."\t".$b[2];
if (($b[3] > 0)&&($b[-1] >= 1000)){
$background{$te_posi}=$b[3];
}}

my $infile1=$name[$index]."_".$p."_25bp_avgCov_unique_mapped";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $locat=$b[0]."\t".$b[1];
$hmd{$locat}=$b[2];
}

my $median;my $k=0;
##focal alternative for distance
my $infile2 = $strain."_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000";

####### control windows for distance
#my $infile2="/dfs7/grylee/yuhenh3/recombination_rate/control_".$strain."_alter_".$alter."_nearest_second_distanc_exclude_TE_overlapped_average_".$control_win_size.".txt";

open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		#my @b = split(" ", $count);
		my @b = split("\t", $count);
#my $chr=$b[0]; my $star=($b[1]+$b[2])/2; my $end=$star; #for between strain comparisons focal strains
###my $chr=$b[5]; my $star=($b[6]+$b[7])/2; my $end=$star; #for between strain comparisons alternative strains
#my $chr=$b[7];my $star=($b[8]+$b[9])/2; my $end=($b[8]+$b[9])/2;#for alternative window
my $chr=$b[-13];my $star=$b[-12]; my $end=$b[-11];

my $te_posi=$chr."\t".$star."\t".$end;
	my @l_enrich; my @r_enrich;my $win; my $dis; my $dis_raw;my @TE_wi_hmd=();
	if (exists $background{$te_posi}){
my $mean_control=$background{$te_posi};

	foreach my $key (keys %hmd) {
my @posi = split("\t", $key);
if ($posi[0] eq $chr){
$win=$posi[1];
if (($win+$hmd_unit >= $star - $win_leng)&&($win <= $end+$win_leng)){
if ($win+$hmd_unit < $star){
push @l_enrich, $hmd{$key}/$mean_control;
}elsif ($win > $end){
push @r_enrich, $hmd{$key}/$mean_control;
}else{
push @TE_wi_hmd, $hmd{$key}/$mean_control;
}
}}
}
my $left_enrich;my $right_enrich;my $avg_enrich;
if ((@l_enrich>=$num_measure)&&(@r_enrich>=$num_measure)){
$left_enrich=mean(@l_enrich);
$right_enrich=mean(@r_enrich);

$avg_enrich=($left_enrich+$right_enrich)/2;
}elsif(@l_enrich>=$num_measure){
	$left_enrich=mean(@l_enrich);
	$right_enrich="NA";
	$avg_enrich=$left_enrich;
}elsif(@r_enrich>=$num_measure){
	$left_enrich="NA";
	$right_enrich=mean(@r_enrich);
	$avg_enrich=$right_enrich;
}else{
	$left_enrich="NA";
	$right_enrich="NA";
	$avg_enrich="NA";
}

if (@TE_wi_hmd>0){
my $avg_TE_hmd=mean(@TE_wi_hmd);
my $median_TE_hmd=median(@TE_wi_hmd);

print OUT "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\t","$left_enrich\t","$right_enrich\t","$avg_enrich\t","$avg_TE_hmd\t","$median_TE_hmd\n";

}else{

print OUT "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\t","$left_enrich\t","$right_enrich\t","$avg_enrich\t","NA\t","NA\n";
}
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
