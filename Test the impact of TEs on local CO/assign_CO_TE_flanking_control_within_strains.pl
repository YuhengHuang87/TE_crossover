#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 5000; my $control_win_size=$win_size*2;

my @strain = ("A6","A7");
foreach my $strain (@strain){
my $outfile1=$strain."_crossover_num_TE_included_window_exclude_TE_overlapped_average_5_10kb_".$win_size.".txt";
open(OUT1, ">$outfile1");

my $outfile2=$strain."_crossover_num_control_included_window_exclude_TE_overlapped_average_".$control_win_size.".txt";
open(OUT2, ">$outfile2");

my $file = '';
my @cross;
my $n=0;my $hit=0;my $loca;
$file =$strain."_cross_output.txt";
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
    #print OUT1 "$a[0]\t","$loca\n";
  }

my $crossover_left=0; my $crossover_location_left='_';
my $crossover_right=0; my $crossover_location_right='_';

my $overlap=0;
my @TE_win_range=();
$file = $strain."_real_pool_coverage_SNP_window_TE_euchromatic_include_ambiguous_".$win_size;
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
        if (@b==9){
    my $chr=$b[1]; my $left=$b[2]-$win_size;my $right=$b[3]+$win_size;
    my $left_bound=$left-$win_size; my $right_bound=$right+$win_size;
    my $range = $chr."\t".$left_bound."\t".$right_bound;
    push @TE_win_range, $range;
}}

  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
            if (@b==9){
    my $chr=$b[1]; my $left=$b[2]-$win_size;my $right=$b[3]+$win_size;
    my $left_bound=$left-$win_size; my $right_bound=$right+$win_size;
    $overlap=0;
    my $focal = $chr."\t".$left_bound."\t".$right_bound;
    for (my $i=0; $i<@TE_win_range; $i++){
      if ($focal ne $TE_win_range[$i]){
      my @d = split("\t", $TE_win_range[$i]);
      if ($d[0] eq $chr){
      if (($d[2]<=$left_bound)||($d[1]>=$right_bound)){}else{
       $overlap++;
      }
      }
    }}
 $crossover_left=0;  $crossover_location_left='_';
 $crossover_right=0;  $crossover_location_right='_';
    if ($overlap==0){
      for (my $j=0; $j<@cross; $j++){
        my @a= split("\t", $cross[$j]);
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
}}


$file = $strain."_real_pool_SNP_coverage_window_control_window_include_ambiguous_".$control_win_size;
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    my $chr=$b[0]; my $left_bound=$b[1];my $right_bound=$b[2];
   my $left=$left_bound+$win_size; my $right=$right_bound-$win_size;
   $crossover_left=0;  $crossover_location_left='_';
$crossover_right=0;  $crossover_location_right='_';

      for (my $j=0; $j<@cross; $j++){
        my @a= split("\t", $cross[$j]);
        if ($chr eq $a[0]){
          if (($a[1] >= $left_bound)&&($a[2] < $left)){
          $crossover_left++;
          $crossover_location_left=$crossover_location_left.$a[1].":".$a[2]."_";
        }elsif (($a[1] > $right)&&($a[2] <= $right_bound)){
          $crossover_right++;
          $crossover_location_right=$crossover_location_right.$a[1].":".$a[2]."_";
      }
    }}
    print OUT2 "$count\t","$crossover_left\t","$crossover_location_left\t","$crossover_right\t","$crossover_location_right\n";
}
}
