#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $strain="A7"; my $alter="A6";
#my $outfile=$strain."_nearest_CO_distance_coverage_SNP_window_TE_euchromatic_include_ambiguous_";
my $outfile="control_nearest_CO_distance_".$strain."_focal_window_depth_TE_include_ambiguous_";

open(OUT, ">$outfile");

my %win_bound;my $window_left;my $chro_window;
my %TE_site;
my $file = ''; my $file1 = ''; 
#$file = "/dfs7/grylee/yuhenh3/recombination_rate/".$strain."_nearest_CO_distance_TE_included_window_exclude_TE_overlapped_average.txt";

$file = "/dfs7/grylee/yuhenh3/recombination_rate/".$strain."_nearest_CO_distance_control_included_window_exclude_TE_overlapped_average.txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    #my $chr= $a[1]; my $left=$a[2]; my $right=$a[3]; #for TE windows
    my $chr= $a[0]; my $left=($a[1]+$a[2])/2; my $right=($a[1]+$a[2])/2;#for control windows
    
    my $TE_call =  $chr."\t".$left."\t".$right;
    my $left_bound = $left+$a[-3]+1;
    my $right_bound = $right+$a[-2]-1;

    my $loca_left = $TE_call."\t"."left"."\t".$left_bound."\t".$left;
    my $loca_right =$TE_call."\t"."right"."\t".$right."\t".$right_bound;
    $TE_site{$loca_left}="";
    $TE_site{$loca_right}="";
  }

my @site_depth=();
my $sum_cover=0;my $loca=''; my $win_loca=''; my $TE_loca='';
my $sum_site_dept=0;

#$file1 = "/dfs7/grylee/yuhenh3/recombination_rate/A6_A6Pool_batch1_2_3_depth_q20_snp_pass.txt";#change for the strains
$file1 = "/dfs7/grylee/yuhenh3/recombination_rate/A7_A7Pool_batch2_3_depth_q20_snp_pass.txt";
  open F, "<$file1" or print "Can't open /$file1\n";
	while(my $count = <F>){
  	chomp($count);
    my @a = split(" ", $count);

    #$sum_site_dept=$a[2]+$a[3]+$a[4];#change for the strains
    $sum_site_dept=$a[2]+$a[3];

    foreach my $key (keys %TE_site) {
      my @b = split("\t",$key);
      if ($a[0] eq $b[0]){
        if (($a[1] > $b[-2])&&($a[1] < $b[-1])){
          $TE_site{$key}=$TE_site{$key}."_".$sum_site_dept;
        }
      }}
  }



  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    #my $chr= $a[1]; my $left=$a[2]; my $right=$a[3];
    my $chr= $a[0]; my $left=int(($a[1]+$a[2])/2); my $right=int(($a[1]+$a[2])/2);#for control windows
    
    my $TE_call =  $chr."\t".$left."\t".$right;
    my $left_bound = $left+$a[-3]+1;
    my $right_bound = $right+$a[-2]-1;
    my $loca_left = $TE_call."\t"."left"."\t".$left_bound."\t".$left;
    my $loca_right =$TE_call."\t"."right"."\t".$right."\t".$right_bound;

    if ((exists $TE_site{$loca_left})&&(exists $TE_site{$loca_right})){
      print "$TE_site{$loca_left}\n";
        my @site_depth_left = split("_",$TE_site{$loca_left});
        splice @site_depth_left, 0, 1;
        my @site_depth_right = split("_",$TE_site{$loca_right});
        splice @site_depth_right, 0, 1;
      if ((@site_depth_left > 0)&&(@site_depth_right > 0)){
        my $mean_depth_left=mean(@site_depth_left);
        my $num_left=scalar(@site_depth_left);

        my $mean_depth_right=mean(@site_depth_right);
        my $num_right=scalar(@site_depth_right);

        print OUT "$count\t","$mean_depth_left\t","$num_left\t","$mean_depth_right\t","$num_right\n";
        }
     }}


my $outfile1=$strain."_nearest_second_distance_TE_included_window_exclude_TE_overlapped_average_".$win_size.".txt"; #windows for TEs
open(OUT1, ">$outfile1");
my $outfile2="control_".$strain."_nearest_second_distanc_exclude_TE_overlapped_average_".$control_win_size.".txt";#control windows
open(OUT2, ">$outfile2");

my $file = '';

my $cov=0; my $snp=0; my $dis=0;
my %TE_win;
$file = "/dfs7/grylee/yuhenh3/recombination_rate/".$strain."_nearest_CO_distance_mean_Mass_TE_HMD_25_1000_10_measures_include_ambiguous.txt";
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    my $family=$b[0]; my $chr=$b[1]; my $left=$b[2];my $right=$b[3]; my $TE_id=$family."\t".$chr."\t".$left."\t".$right;
    my @left_right_dis=(abs($b[8]),abs($b[9]),abs($b[10]),abs($b[11]));
    my @sorted_distance = sort { $a <=> $b } @left_right_dis;
    if ($sorted_distance[0]==abs($b[9])){
        $cov=$b[12];$snp=$b[13];
    }elsif($sorted_distance[0]==abs($b[10])){
        $cov=$b[14];$snp=$b[15];
    }
    $TE_win{$TE_id}=$b[4]."\t".$sorted_distance[0]."\t".$sorted_distance[1]."\t".$cov."\t".$snp."\t".$b[-6]."\t".$b[-5]."\t".$b[-4]."\t".$b[-3]."\t".$b[-2]."\t".$b[-1];
        print OUT1 "$TE_id\t","$TE_win{$TE_id}\n";
  }

  $file = "control_nearest_CO_distance_".$strain."_focal_window_depth_TE_include_ambiguous_";
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    my $chr=$b[0]; my $left=$b[1];my $right=$b[2]; my $win_id=$chr."\t".$left."\t".$right;

    my @left_right_dis=(abs($b[3]),abs($b[4]),abs($b[5]),abs($b[6]));
    my @sorted_distance = sort { $a <=> $b } @left_right_dis;
    if ($sorted_distance[0]==abs($b[4])){
        $cov=$b[-4];$snp=$b[-3];
    }elsif($sorted_distance[0]==abs($b[5])){
        $cov=$b[-2];$snp=$b[-1];
    }
    print OUT2 "$win_id\t","$sorted_distance[0]\t","$sorted_distance[1]\t","$cov\t","$snp\n";
  }

sub mean {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}
