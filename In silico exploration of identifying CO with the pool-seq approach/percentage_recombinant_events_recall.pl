#!/usr/bin/perl
use strict;
use warnings;

my $total_CO=2000;
my @depth=("0.1","0.25","0.5","1","2","3","4","6","8");
#my @depth=("3");
my @chromo=("CM010569.1","CM010570.1","CM010571.1","CM010572.1","CM010573.1");
my @cutoff=("0","1500","2000","2500");
    	foreach my $cf (@cutoff){
my $outfile=$cf."_depth_event_identified_read_num.txt";
open(OUT, ">$outfile");
    	foreach my $depth (@depth){

my $outfile1=$depth."_event_identified_".$cf;
open(OUT1, ">$outfile1");

my $outfile2=$depth."_event_unidentified_".$cf;
open(OUT2, ">$outfile2");
my $k=0; my $sum_ident=0;

my %share;my %event_total; my %event_found;
my $infile2=$depth."_co_gc_shared_each_read_fragment_cutoff_".$cf."_10_10_10_2";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
my @b = split("\t", $count);
    my @a = split("_", $b[0]);
    my @d = split("S", $a[0]);
my $id = $d[1]."\t".$a[-1]."\t".$b[0];
	$share{$id}=$b[1]."\t".$b[-4]."\t".$b[-3];
	$k++;
}

foreach my $chromo (@chromo){
 my $i=0;my $j=0;my $m=0;
my $infile1="A6_A4_recombine_individual_arms_".$chromo."_location";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
        $m=0;
        foreach my $key (keys %share){
            my @a = split("\t", $key);
            if (($a[1] eq $chromo)&&($a[0] eq $b[0])){
                my @c = split("\t", $share{$key});
                if (($c[1]<=$b[2])&&($c[2]>=$b[2])){
                    $m++;
                }
            }}

if ($m>0){
$i+=1;
print OUT "$depth\t","$count\t","$m\n";
print OUT1 "$count\n";
}else{
$j++;
print OUT2 "$count\n";
}
		}
close FILE1;

my $sum_event=$i+$j;
$sum_ident=$sum_ident+$i;
        }
my $prop=$sum_ident/$total_CO;
print "$depth\t","$cf\t","$sum_ident\t","$prop\n";
}}