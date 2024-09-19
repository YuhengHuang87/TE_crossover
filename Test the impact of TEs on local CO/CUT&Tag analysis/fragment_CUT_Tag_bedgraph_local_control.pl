#!/usr/bin/perl
use strict;
use warnings;
my $site_num=25;
my @ind=("A6_1","A6_2","A6_4");
#my @ind=("A7_2","A7_3","A7_4");

#my @euch=("CM010569.1_144265_18689904","CM010570.1_251300_19760175","CM010571.1_8030989_23829870","CM010572.1_164091_18449201","CM010573.1_8729782_31194759");
#my @euch=("CM010576.1_204176_18599896","CM010577.1_134798_19487247","CM010578.1_9035199_24726584","CM010579.1_292296_18448828","CM010580.1_9655508_32200062");

foreach my $ind (@ind){
my $outfile=$ind."_".$site_num."bp_avgCov_unique_mapped";
open(OUT, ">$outfile");

my $sum_cover=0; my $sum_site=0; my $chr=''; my $i; my $win_start; my $chr_win;
my $win_cov; my $win_num=0;


my $infile=$ind."_bowtie2.fragments.bedgraph";
open(FILE1,"<$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
my @b=split("\t",$count);
#if (($b[0] eq "CM010541.1")||($b[0] eq "CM010542.1")||($b[0] eq "CM010543.1")||($b[0] eq "CM010544.1")||($b[0] eq "CM010545.1")){
if (($b[0] eq "CM010569.1")||($b[0] eq "CM010570.1")||($b[0] eq "CM010571.1")||($b[0] eq "CM010572.1")||($b[0] eq "CM010573.1")){
#if (($b[0] eq "CM010576.1")||($b[0] eq "CM010577.1")||($b[0] eq "CM010578.1")||($b[0] eq "CM010579.1")||($b[0] eq "CM010580.1")){
my $leng=$b[2]-$b[1];
my $cov=$b[3];
if ($b[0] eq $chr){
for ($i = 1; $i <= $leng; $i++){
$sum_site++;
$sum_cover+=$cov;
if ($sum_site == $site_num){
$win_cov=$sum_cover/$sum_site;
print OUT "$chr_win\t","$win_cov\n";
$win_num++;
$sum_cover=0;
$sum_site=0;
$win_start=$b[1]+$i;
$chr_win=$b[0]."\t".$win_start;
}}
}else{
if($chr eq ''){
$win_start=$b[1];
$chr_win=$b[0]."\t".$win_start;
$sum_cover=0;
$sum_site=0;
$chr = $b[0];
}else{
if ($sum_site>0){
$win_cov=$sum_cover/$sum_site;
print OUT "$chr_win\t","$win_cov\n";
}
$win_num++;
$win_start=$b[1];
$chr_win=$b[0]."\t".$win_start;
$sum_cover=0;
$sum_site=0;
$chr = $b[0];
}
for ($i = 1; $i <= $leng; $i++){
$sum_site++;
$sum_cover+=$cov;
if ($sum_site == $site_num){
$win_cov=$sum_cover/$sum_site;
print OUT "$chr_win\t","$win_cov\n";
$win_num++;
$sum_cover=0;
$sum_site=0;
$win_start=$b[1]+$i;
$chr_win=$b[0]."\t".$win_start;
}}
}}
}}

my @name=("A6","A7","A4");
my @combine=('1_2_4','2_3_4');
my $strain="A6"; my $alter="A7";
my $index=1; #need to match strain
my $rep=0;
my $win_size= 5000; my $dis1=30000;my $dis2=40000;
my @index_pair1=split("_", $combine[$index]);
my $p= $index_pair1[$rep];


my $i=0;my $j=0;my $k=0;my %TE_ref;

my $outfile1="HMD_unique_mapped_local_nearest_distance_alter_".$alter."_".$p.".txt";
open(OUT, ">$outfile1");

my $infile1=$strain."_crossover_nearest_distance_coverage_SNP_TE_left_right_window_mass_5000_exclude_distance_5000";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $chr=$b[-13];my $left=$b[-12]; my $right=$b[-11];

#my $chr=$b[0];my $left=$b[1]; my $right=$b[2];
my $locat=$chr."\t".$left."\t".$right;
$i++;
$TE_ref{$locat}=$b[0];
}
#}}

foreach my $key (keys %TE_ref) {
my @posi = split("\t", $key);
my @local_hmd;my $TE_num=0;
my $te_l; my $te_r;
my $infile2=$name[$index]."_".$p."_25bp_avgCov_unique_mapped";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
	$j++;
if ($b[0] eq $posi[0]){
$te_l=$posi[1];$te_r=$posi[2]; my $win=$b[1];
if ((($win > $te_l - 40000)&&($win+25 < $te_l-20000))||(($win+25 < $te_r + 40000)&&($win > $te_r+20000))){
push @local_hmd, $b[2];
	$TE_num++;
}
}}

my $hmd_median=median(@local_hmd);
my $hmd_mean=mean(@local_hmd);
print OUT "$key\t","$hmd_mean\t","$hmd_median\t","$TE_num\n";
}
#print "$name\t","$i\t","$j\n";


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