#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my $outfile1="aln_A6_benchmark_pool_pbmm2_sorted_que20_4kb.sam";
open(OUT1, ">$outfile1");

my $file = '';
my %read; my $n=0;my $hit=0;

$file = "aln_A6_benchmark_pool_sorted.sam";
#$file = "aln_A4_benchmark_pool_sorted.sam";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    if($count =~ m/@/){
      print OUT1 "$count\n";
    }else{
    my @a = split("\t", $count);
    my $chr=$a[2]; my $left=$a[3]; my $right=$a[3]+2000;
    #if ((($chr eq "CM010541.1")&&($left>206715)&&($right<18899017))||(($chr eq "CM010542.1")&&($left>123450)&&($right<19558299))||(($chr eq "CM010543.1")&&($left>8597017)&&($right<24417585))||(($chr eq "CM010544.1")&&($left>173141)&&($right<18439213))||(($chr eq "CM010545.1")&&($left>10069791)&&($right<32445046))){ #Dmel for A4
    if ((($chr eq "CM010569.1")&&($left>144265)&&($right<18689904))||(($chr eq "CM010570.1")&&($left>251300)&&($right<19760175))||(($chr eq "CM010571.1")&&($left>8030989)&&($right<23829870))||(($chr eq "CM010572.1")&&($left>164091)&&($right<18449201))||(($chr eq "CM010573.1")&&($left>8729782)&&($right<31194759))){ #Dmel for A6
    #if ((($chr eq "CM010576.1")&&($left>204176)&&($right<18599896))||(($chr eq "CM010577.1")&&($left>134798)&&($right<19487247))||(($chr eq "CM010578.1")&&($left>9035199)&&($right<24726584))||(($chr eq "CM010579.1")&&($left>292296)&&($right<18448828))||(($chr eq "CM010580.1")&&($left>9655508)&&($right<32200062))){ #Dmel for A7
    if ($a[4]>=20){
      my @d=split("",$a[9]);
      if (@d>=4000){
          print OUT1 "$count\n";
          $n++;
    }
  }
}}
}
print "$hit\t","$n\n";
