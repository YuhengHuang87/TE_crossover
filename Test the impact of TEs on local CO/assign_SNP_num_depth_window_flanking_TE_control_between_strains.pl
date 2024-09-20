#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 5000; my $dis = 5000;
my $strain="A7"; my $alter="A6";
my $outfile="TE_focal_".$strain."_alter_".$alter."_SNP_coverage_window_TE_exclude_".$win_size."_exclude_distance_".$dis;
open(OUT, ">$outfile");

my %win_bound;my $window_left;my $chro_window;

my $file = ''; my %TE_site;
$file = $strain."_focal_".$alter."_location_include_nearby_ambig_clean_TE_".$win_size."_exclude_distance_".$dis;
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    my $chr= $a[-3]; my $left=$a[-2]; my $right=$a[-1];
    my $TE_call =  $count;
   my $left_bound = $left-$win_size+1;
    my $right_bound = $right+$win_size-1;
  
    for (my $j=$left_bound; $j<=$left;$j++){
    my $loca = $chr."\t".$j;
    $TE_site{$loca}=$TE_call.":"."left";
  }
    for (my $j=$right; $j<=$right_bound;$j++){
    my $loca = $chr."\t".$j;
    $TE_site{$loca}=$TE_call.":"."right";
  }
}

my @site_depth=();
my $sum_cover=0;my $loca=''; my $win_loca=''; my $TE_loca='';
my $sum_site_dept=0;

$file = "A6_A6Pool_batch1_2_3_depth_q20_snp_pass.txt";
#$file = "A7_A7Pool_batch2_3_depth_q20_snp_pass.txt";

  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split(" ", $count);
    $loca = $a[0]."\t".$a[1];

    $sum_site_dept=$a[2]+$a[3]+$a[4];
    if (exists $TE_site{$loca}){
              my @b = split(":",$win_loca);
        if ($TE_site{$loca} ne $win_loca){
      if (@site_depth>0){
        my $mean_depth=mean(@site_depth);
        my $num=scalar(@site_depth);
        if ($b[0] ne $TE_loca){
        print OUT "\n";
        print OUT "$b[0]\t","$mean_depth\t","$num\t";
        }elsif($b[0] eq $TE_loca){
        print OUT "$mean_depth\t","$num";
        }
      }
    $TE_loca=$b[0];
    @site_depth=();
    $win_loca=$TE_site{$loca};
    push @site_depth,$sum_site_dept;

    }else{
    push @site_depth,$sum_site_dept;

    }
    }
    }
    
  if (@site_depth>0){
        my $mean_depth=mean(@site_depth);
        my $num=scalar(@site_depth);
        print OUT "$mean_depth\t","$num\n";
      }

### depending on which strain as focal, swap the input files 1 and 2 
my $file2 = "/dfs7/grylee/yuhenh3/recombination_rate/A6_A6Pool_batch1_2_3_depth_q20_snp_pass.txt";
my $file1 = "/dfs7/grylee/yuhenh3/recombination_rate/A7_A7Pool_batch2_3_depth_q20_snp_pass.txt";

my $outfile1="control_focal_".$strain."_alter_".$alter."_SNP_coverage_window_TE_exclude_".$dis.".txt"; 
open(OUT1, ">$outfile1");

my $file = '';
my @win_control_focal; my @win_control_alter; 
$file = $strain."_focal_".$alter."_location_include_nearby_ambig_clean_control_".$dis.".txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    my $chr = $a[-3]; my $left=$a[-2]; my $right=$a[-1];
    my $loca = $chr."\t".$left."\t".$right;
    push @win_control_alter, $loca;

    $chr = $a[1]; $left=$a[2]; $right=$a[3];
    $loca = $chr."\t".$left."\t".$right;
    push @win_control_focal, $loca;
  }

for (my $k=0; $k<@win_control_focal; $k+=1){
my @b = split("\t", $win_control_focal[$k]);
my $chro=$b[0];my $window_left=$b[1];my $window_right = $b[2];
@site_depth=();

  open I, "<$file1" or print "Can't open /$file\n";#file1
	while(my $count = <I>){
  	chomp($count);
    my @a = split(" ", $count);
    $loca = $a[0]."\t".$a[1];
    $sum_site_dept=$a[2]+$a[3]+$a[4]; 
    if ($a[0] eq $chro){
    if (($a[1]>=$window_left)&&($a[1]<=$window_right)){
      push @site_depth,$sum_site_dept;
    }}}
      if (@site_depth>0){
        my $mean_depth=mean(@site_depth);
        my $num=scalar(@site_depth);
        print OUT1 "$win_control_focal[$k]\t","$mean_depth\t","$num\t";
      }
      
      @site_depth=();
      @b = split("\t", $win_control_alter[$k]);
      $chro=$b[0];$window_left=$b[1];$window_right = $b[2];

    open I, "<$file2" or print "Can't open /$file\n";#file2
	while(my $count = <I>){
  	chomp($count);
    my @a = split(" ", $count);
    $loca = $a[0]."\t".$a[1];
    $sum_site_dept=$a[2]+$a[3]+$a[4]; 
    if ($a[0] eq $chro){
    if (($a[1]>=$window_left)&&($a[1]<=$window_right)){
      push @site_depth,$sum_site_dept;
    }}}
      if (@site_depth>0){
        my $mean_depth=mean(@site_depth);
        my $num=scalar(@site_depth);
        print OUT1 "$win_control_alter[$k]\t","$mean_depth\t","$num\n";
      }else{
        print OUT1 "\n";
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
