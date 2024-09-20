#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 5000; my $dis=5000;
my $strain="A7";my $alter="A6";
my $outfile=$strain."_focal_".$alter."_location_include_nearby_ambig_clean_TE_".$win_size."_exclude_distance_".$dis;
open(OUT, ">$outfile");

my $outfile1=$strain."_focal_".$alter."_location_include_nearby_ambig_clean_control_".$dis.".txt"; 
open(OUT1, ">$outfile1");



my $file = ''; my %focal_posi;
my @focal_te; my @alter_te;
my $n=0;my $hit=0;
$file = $alter."focal_5000bp_distance_".$strain."_location_TE_".$win_size.".txt";
#$file = $strain_focal."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic.txt";

  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    my $chr = $a[-3]; my $left=$a[-2]; my $right=$a[-1];
    my $focal_loca = $chr."\t".$left."\t".$right;
          push @focal_te,$focal_loca;
    my $alter_loca = $a[1]."\t".$a[2]."\t".$a[3];
      push @alter_te,$alter_loca;
    }

$file = $strain."focal_5000bp_distance_".$alter."_location_TE_".$win_size.".txt";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
  chomp($count);
  my @b = split("\t", $count);
  my $alter_chr= $b[-3]; my $alter_left=$b[-2]; my $alter_right=$b[-1];
  my $focal_chr= $b[1]; my $focal_left=$b[2]; my $focal_right=$b[3];
        my $overlap=0;
            for (my $j=0; $j<@focal_te; $j++){
        my @a= split("\t", $focal_te[$j]);
        my $te_left=$a[1]-$dis; my $te_right = $a[2]+$dis;
        if ($a[0] eq $focal_chr){ 
        if (($te_left >= $focal_right)||($te_right <= $focal_left)){
        }else{
          $overlap++;
        }
      }}
            for (my $j=0; $j<@alter_te; $j++){
        my @a= split("\t", $alter_te[$j]);
        my $te_left=$a[1]-$dis; my $te_right = $a[2]+$dis;
        if ($a[0] eq $alter_chr){
        if (($te_left >= $alter_right)||($te_right <= $alter_left)){
        }else{
          $overlap++;
        }
      }}
      if ($overlap==0){
    print OUT "$b[0]\t","$b[1]\t","$b[2]\t","$b[3]\t","$b[4]\t","$b[5]\t","$b[6]\n";
    $n++;
  }
}
print "$n\n";


$file = $strain."_focal_control_alter_".$alter."_location_window_10000.txt";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
  chomp($count);
  my @b = split("\t", $count);
  my $alter_chr= $b[-3]; my $alter_left=$b[-2]; my $alter_right=$b[-1];
  #my $focal_chr= $b[1]; my $focal_left=$b[2]; my $focal_right=$b[3];
        my $overlap=0;
        for (my $j=0; $j<@alter_te; $j++){
        my @a= split("\t", $alter_te[$j]);
        my $te_left=$a[1]-$dis; my $te_right = $a[2]+$dis;
        if ($a[0] eq $alter_chr){
        if (($te_left > $alter_right)||($te_right < $alter_left)){
        }else{
          $overlap++;
        }
      }}
if ($overlap==0){
    print OUT1 "$count\n";
  }
}


