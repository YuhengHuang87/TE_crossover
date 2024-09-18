#!/usr/bin/perl
use strict;
use warnings;

my @locu_ref;
my @rec_ref;

my $a=0; my $t=0;
my $strain="A6";#or A7
my $outfile="recombination_rate_".$strain.".txt";
open(OUT, ">$outfile");

my $infile1="Comeron_2012_PLoSGent_inR6.txt";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $left = $b[1]-50000; my $right = $b[1]+50000;
my $locat=$b[0]."\t".$left."\t".$right;
push @locu_ref,$locat;
push @rec_ref,$b[2];
$t++;
}

for (my $k=0; $k<@locu_ref; $k++){
	print OUT "$locu_ref[$k]\t","$rec_ref[$k]\t";
my @posi = split("\t", $locu_ref[$k]);
my $chr=$posi[0];
my $i=0; my $j=0;
my $ref_str=''; my $ref_end=''; my $que_str=''; my $que_end='';
my $infile2="dmel_6.43_refer_A7_query_asm10_update.var.txt";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
if($b[0] eq "V"){
if ($b[1] eq $chr){
$ref_end=$b[2]; #need to change
$que_end=$b[9];

#$ref_end=$b[9];
#$que_end=$b[2];

my $te_q_s; my $te_q_e;
if ($ref_str ne ""){
if(($posi[1]>=$ref_str)&&($posi[1]<=$ref_end)){
my $up_dis_s=$posi[1]-$ref_str;
$te_q_s=$que_str+$up_dis_s;
$i++;
print OUT "$b[-4]\t","$te_q_s\t";
}
if(($posi[2]>=$ref_str)&&($posi[2]<=$ref_end)){
my $up_dis_e=$posi[2]-$ref_str;
$te_q_e=$que_str+$up_dis_e;
$j++;
print OUT "$te_q_e\t";
}
if (($i==1)&&($j==1)){
last;
}
}
$ref_str=$b[3];#need to change
$que_str=$b[10];

#$ref_str=$b[10];
#$que_str=$b[3];
}}}
print OUT "\n";
}

my $win_size=50000;
my $outfile1=$strain."_Dm6_rates_crossover_TE_num_control_window_weight_".$win_size.".txt";
open(OUT1, ">$outfile1");

my $file = '';
my @cross;
my $n=0;my $hit=0;my $loca;
$file = $strain."_cross_output.txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @a = split("\t", $count);
    $loca = $a[0]."\t".$a[1]."\t".$a[2];
      push @cross,$loca;
  }

my $crossover=0; my $peak_area=0;
my $window_left;my $window_right;
$file = "recombination_rate_".$strain.".txt";
  open I, "<$file" or print "Can't open /$file\n";
  while(my $count = <I>){
    chomp($count);
    my @b = split("\t", $count);
    $crossover=0;
    #my $chr=$b[0]; my $left=$b[1];my $right=$b[2];
    my $chr=$b[-2]; my $left=$b[-1]-$win_size;my $right=$b[-1]+$win_size;

    if ((($chr eq "CM010569.1")&&($left>144265)&&($right<18689904))||(($chr eq "CM010570.1")&&($left>251300)&&($right<19760175))||(($chr eq "CM010571.1")&&($left>8030989)&&($right<23829870))||(($chr eq "CM010572.1")&&($left>164091)&&($right<18449201))||(($chr eq "CM010573.1")&&($left>8729782)&&($right<31194759))){ #Dmel for A6
    #if ((($chr eq "CM010576.1")&&($left>204176)&&($right<18599896))||(($chr eq "CM010577.1")&&($left>134798)&&($right<19487247))||(($chr eq "CM010578.1")&&($left>9035199)&&($right<24726584))||(($chr eq "CM010579.1")&&($left>292296)&&($right<18448828))||(($chr eq "CM010580.1")&&($left>9655508)&&($right<32200062))){ #Dmel for A7
      for (my $j=0; $j<@cross; $j++){
        my @a= split("\t", $cross[$j]);
        if ($chr eq $a[0]){
        if (($a[1] >= $left)&&($a[2] <= $right)){
          $crossover++;
        }elsif (($a[1] < $left)&&($a[2] > $left)){
          $crossover=$crossover+0.5;
        }elsif(($a[1] < $right)&&($a[2] > $right)){
          $crossover=$crossover+0.5;
        }
      }
    }
    print OUT1 "$count\t","$crossover\n";
  }
}

my $outfile2="CO_broad_".$strain."_focal_window_depth.txt";
open(OUT2, ">$outfile2");

my %win_bound;my $window_left;my $chro_window;
my %win_site;
my $file = ''; my $file1 = ''; 

$file = $strain."_Dm6_rates_crossover_TE_num_control_window_weight_".$win_size.".txt";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @b = split("\t", $count);
    my $chr= $b[3]; my $left=$b[4]-$win_size;my $right=$b[4]+$win_size;
    my $win_call =  $chr."\t".$left."\t".$right;
    $win_site{$win_call}="";
  }

my @site_depth=();
my $sum_cover=0;my $loca=''; my $win_loca=''; my $TE_loca='';
my $sum_site_dept=0;

$file1 = $strain."_assembly_".$strain."_pool_read_pass_snp_".$snp_num"._depth.txt";
  open F, "<$file1" or print "Can't open /$file1\n";
	while(my $count = <F>){
  	chomp($count);
    my @a = split(" ", $count);
    #$sum_site_dept=$a[2]+$a[3]+$a[4];#change for the strains
    $sum_site_dept=$a[2]+$a[3];

    foreach my $key (keys %win_site) {
      my @b = split("\t",$key);
      if ($a[0] eq $b[0]){
        if (($a[1] > $b[-2])&&($a[1] < $b[-1])){
          $win_site{$key}=$win_site{$key}."_".$sum_site_dept;
        }
      }}
  }

  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    my @b = split("\t", $count);
    my $chr= $b[3]; my $left=$b[4]-$win_size;my $right=$b[4]+$win_size;
    my $win_call =  $chr."\t".$left."\t".$right;
   
    if (exists $win_site{$win_call}){
      print "$win_site{$win_call}\n";
        my @site_depth = split("_",$win_site{$win_call});
        splice @site_depth, 0, 1;
        
      if (@site_depth > 0){
        my $mean_depth=mean(@site_depth)/2;
        my $num=scalar(@site_depth);
        my $percent_CO = $b[5]/$mean_depth;
        print OUT2 "$count\t","$percent_CO\t","$mean_depth\t","$num\n";
        }
     }}


sub mean {
my @array = @_; # save the array passed to this function
my $sum; # create a variable to hold the sum of the array's values
foreach (@array) { $sum += $_; } # add each element of the array
# to the sum
return $sum/@array; # divide sum by the number of elements in the
# array to find the mean
}