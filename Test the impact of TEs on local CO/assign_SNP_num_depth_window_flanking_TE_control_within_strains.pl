#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $win_size= 50000;
my $strain="A6"; 
my $outfile=$strain."_real_pool_coverage_SNP_window_TE_euchromatic_include_ambiguous_".$win_size;
open(OUT, ">$outfile");

my %win_bound;my $window_left;my $chro_window;

my $file = ''; my %TE_site;
$file = $strain."_exclu_within_5000bp_distance_nearbyall_TE_euchromatic.txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    my $chr= $a[1]; my $left=$a[2]; my $right=$a[3];
    
    my $TE_call =  $a[0]."\t".$chr."\t".$left."\t".$right."\t".$a[4];
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

$file = "A6_A6Pool_batch1_2_3_depth_q20_snp_pass.txt";#change the strains for different files containing the depth info

  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split(" ", $count);
    $loca = $a[0]."\t".$a[1];
        #$sum_site_dept=$a[2];

    $sum_site_dept=$a[2]+$a[3]+$a[4];#change for the strains
    #$sum_site_dept=$a[2]+$a[3];


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


my $control_size=$win_size*2;
my $outfile1=$strain."_real_pool_SNP_coverage_window_control_window_include_ambiguous_".$control_size;
open(OUT1, ">$outfile1");

my @euch=("CM010569.1_144265_18689904","CM010570.1_251300_19760175","CM010571.1_8030989_23829870","CM010572.1_164091_18449201","CM010573.1_8729782_31194759");
#my @euch=("CM010576.1_204176_18599896","CM010577.1_134798_19487247","CM010578.1_9035199_24726584","CM010579.1_292296_18448828","CM010580.1_9655508_32200062");

my @site_depth_1=();my @site_depth_2=();
my $sum_cover=0;my $loca=''; my $win_loca=''; my $overlap=0;

my $sum_site_dept=0;
my @win_bound;my $window_left;my $chro_window;
foreach my $euch (@euch){
my @b = split("_", $euch);
my $chro=$b[0];
my $start=$b[1];
my $end = $b[2];
for (my $i=0; $i<25000000; $i+=$win_size){
$overlap=0;
  $window_left=$start+$i;
  $window_right=$window_left+$win_size;
  if ($window_right<=$end){
  $chro_window=$chro."\t".$window_left."\t".$window_right;
for (my $j=0; $j<@TE_site; $j+=1){
    my @TE = split("\t", $TE_site[$j]);
    if ($TE[0] eq $chro){
        if (($TE[1]-$dis>$window_right)||($TE[2]+$dis<$window_left)){}else{
            $overlap++;
        }}}
if ($overlap==0){
  push @win_bound, $chro_window;
}
  }}}

for (my $k=0; $k<@win_bound; $k+=1){
    #print "$win_bound[$k]\n";
my @b = split("\t", $win_bound[$k]);
my $chro=$b[0];
my $win_1_left=$b[1]; my $win_1_right = $b[1]+$win_size/2;
my $win_2_left=$b[2]-$win_size/2; my $win_2_right = $b[2];

  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split(" ", $count);
    $loca = $a[0]."\t".$a[1];
    $sum_site_dept=$a[2]+$a[3]+$a[4];#change it for A7
    #$sum_site_dept=$a[2]+$a[3];

    if ($a[0] eq $chro){
    if (($a[1]>=$win_1_left)&&($a[1]<=$win_1_right)){
      push @site_depth_1,$sum_site_dept;
    }
  if (($a[1]>=$win_2_left)&&($a[1]<=$win_2_right)){
      push @site_depth_2,$sum_site_dept;
    }
  }}

      if ((@site_depth_1>0)&&(@site_depth_2>0)){
        my $mean_depth_1=mean(@site_depth_1);my $num_1=scalar(@site_depth_1);
        my $mean_depth_2=mean(@site_depth_2);my $num_2=scalar(@site_depth_2);
        print OUT "$win_bound[$k]\t","$mean_depth_1\t","$num_1\t","$mean_depth_2\t","$num_2\n";
      }
      @site_depth_1=();@site_depth_2=();
    }


sub mean {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}
