#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 5000; my $control_win_size=$win_size*2;

my @strain = ("A6","A7");
foreach my $strain (@strain){
my $outfile1=$strain."_nearest_CO_distance_TE_included_window_exclude_TE_overlapped_average.txt";
open(OUT1, ">$outfile1");

my $outfile2=$strain."_nearest_CO_distance_control_included_window_exclude_TE_overlapped_average.txt";
open(OUT2, ">$outfile2");

my $file = '';
my @cross;
my $n=0;my $hit=0;my $loca;
$file = $strain."_cross_output.txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    if (@a==6){
    $loca = $a[1]."\t".$a[-2]."\t".$a[-1];
  }else{
    $loca = $a[1]."\t".$a[2]."\t".$a[3];
  }
      push @cross,$loca;
  }

my $crossover_left_dis=99999999; my $crossover_location_left='_';
my $crossover_right_dis=99999999; my $crossover_location_right='_';

my $overlap=0;
my @TE_win_range=();
$file =  $strain."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic.txt";
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    my $chr=$b[1]; my $left=$b[2];my $right=$b[3];
    $overlap=0;
 $crossover_left_dis=99999999;  $crossover_location_left='_';
 $crossover_right_dis=99999999;  $crossover_location_right='_';
    my @left_dis=(); my @right_dis=();
      for (my $j=0; $j<@cross; $j++){
        my @a= split("\t", $cross[$j]);
        if ($chr eq $a[0]){
          if ($a[2] < $left){
            my $distance = $a[2]-$left;
            push @left_dis,$distance;
        }elsif ($a[1] > $right){
            my $distance = $a[1]-$right;
            push @right_dis,$distance;
      }
    }}
    if ((@left_dis>1)&&(@right_dis>1)){
    my @sort_left = sort { $b <=> $a } @left_dis;
    my @sort_right = sort { $a <=> $b } @right_dis;
    print OUT1 "$count\t","$sort_left[1]\t","$sort_left[0]\t","$sort_right[0]\t","$sort_right[1]\n";
}}


$file = $strain."_real_pool_SNP_coverage_window_control_window_include_ambiguous_".$control_win_size;
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    my $chr=$b[0]; my $left_bound=$b[1];my $right_bound=$b[2];
   my $mid = int(($left_bound+$right_bound)/2);
 $crossover_left_dis=99999999;  $crossover_location_left='_';
 $crossover_right_dis=99999999;  $crossover_location_right='_';
    my @left_dis=(); my @right_dis=();
      for (my $j=0; $j<@cross; $j++){
        my @a= split("\t", $cross[$j]);
        if ($chr eq $a[0]){
          if ($a[2] < $mid){
            my $distance = $a[2]-$mid;
            push @left_dis,$distance;
        }elsif ($a[1] > $mid){
            my $distance = $a[1]-$mid;
            push @right_dis,$distance;
      }
    }}
        if ((@left_dis>1)&&(@right_dis>1)){
    my @sort_left = sort { $b <=> $a } @left_dis;
    my @sort_right = sort { $a <=> $b } @right_dis;
    print OUT2 "$b[0]\t","$b[1]\t","$b[2]\t","$sort_left[1]\t","$sort_left[0]\t","$sort_right[0]\t","$sort_right[1]\n";
}}
}
