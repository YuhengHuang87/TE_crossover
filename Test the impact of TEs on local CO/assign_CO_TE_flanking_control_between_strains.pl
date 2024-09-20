#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 5000; my $control_win_size=$win_size*2; 
my $strain="A7"; my $alter="A6";

my $outfile1="focal_".$strain."_alter_".$alter."_crossover_num_TE_included_window_exclude_TE_overlapped_average_".$win_size.".txt";
open(OUT1, ">$outfile1");

my $outfile2="focal_".$strain."_alter_".$alter."_crossover_num_control_included_window_exclude_TE_overlapped_average_".$control_win_size.".txt";
open(OUT2, ">$outfile2");

my $file = '';
my @cross_focal; my @cross_alter;
my $loca;
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
      push @cross_focal,$loca;
  }

$file = $alter."_cross_output.txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    if (@a==6){
    $loca = $a[1]."\t".$a[-2]."\t".$a[-1];
  }else{
    $loca = $a[1]."\t".$a[2]."\t".$a[3];
  }
      push @cross_alter,$loca;
  }

my $crossover_left=0; my $crossover_location_left='_';
my $crossover_right=0; my $crossover_location_right='_';

my $overlap=0;
my @TE_win_range=();
$file = "TE_focal_".$strain."_alter_".$alter."_SNP_coverage_window_TE_exclude_".$win_size;
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
        if (@b==11){
    my $chr=$b[4]; my $left=$b[5];my $right=$b[6];
    my $left_bound=$left-$win_size; my $right_bound=$right+$win_size;
    $overlap=0;
    my $focal = $chr."\t".$left_bound."\t".$right_bound;
 $crossover_left=0;  $crossover_location_left='_';
 $crossover_right=0;  $crossover_location_right='_';
      for (my $j=0; $j<@cross_alter; $j++){
        my @a= split("\t", $cross_alter[$j]);
        if ($chr eq $a[0]){
          if (($a[1] > $left_bound)&&($a[2] < $left)){
          $crossover_left++;
          $crossover_location_left=$crossover_location_left.$a[1].":".$a[2]."_";
        }elsif (($a[1] > $right)&&($a[2] < $right_bound)){
          $crossover_right++;
          $crossover_location_right=$crossover_location_right.$a[1].":".$a[2]."_";
      }
    }}
    print OUT1 "$count\t","$crossover_left\t","$crossover_location_left\t","$crossover_right\t","$crossover_location_right\n";
  }else{
    print "$count\n";
  }
}


$file = "control_focal_".$strain."_alter_".$alter."_SNP_coverage_window_TE_exclude_".$control_win_size.".txt";
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    if (@b == 10){
    my $chr=$b[0]; my $left_bound=$b[1];my $right_bound=$b[2];my $left=$left_bound+$win_size; my $right=$right_bound-$win_size;
    $crossover_left=0;  $crossover_location_left='_';$crossover_right=0;  $crossover_location_right='_'; 

      for (my $j=0; $j<@cross_focal; $j++){
        my @a= split("\t", $cross_focal[$j]);
        if ($chr eq $a[0]){
          if (($a[1] >= $left_bound)&&($a[2] < $left)){
          $crossover_left++;
          $crossover_location_left=$crossover_location_left.$a[1].":".$a[2]."_";
        }elsif (($a[1] > $right)&&($a[2] <= $right_bound)){
          $crossover_right++;
          $crossover_location_right=$crossover_location_right.$a[1].":".$a[2]."_";
      }
    }}
print OUT2 "$count\t","$crossover_left\t","$crossover_location_left\t","$crossover_right\t","$crossover_location_right\t";
 $chr=$b[5];  $left_bound=$b[6]; $right_bound=$b[7]; $left=$left_bound+$win_size;  $right=$right_bound-$win_size;
    $crossover_left=0;  $crossover_location_left='_';$crossover_right=0;  $crossover_location_right='_'; 

      for (my $j=0; $j<@cross_alter; $j++){
        my @a= split("\t", $cross_alter[$j]);
        if ($chr eq $a[0]){
          if (($a[1] >= $left_bound)&&($a[2] < $left)){
          $crossover_left++;
          $crossover_location_left=$crossover_location_left.$a[1].":".$a[2]."_";
        }elsif (($a[1] > $right)&&($a[2] <= $right_bound)){
          $crossover_right++;
          $crossover_location_right=$crossover_location_right.$a[1].":".$a[2]."_";
      }
    }}
print OUT2 "$crossover_left\t","$crossover_location_left\t","$crossover_right\t","$crossover_location_right\n";
}
  }

