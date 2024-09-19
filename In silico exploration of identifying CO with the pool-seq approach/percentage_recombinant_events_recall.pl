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

my $event_num=10; my $depth=3;
my @ID=("CM010569.1","CM010570.1","CM010571.1","CM010572.1","CM010573.1");
my $outfile3=$depth."_co_region_false_negative_rate_nonoverlap_".$event_num.".txt";
open(OUT3, ">$outfile3");

foreach my $ID (@ID){
my @win_allele=();my @win_read=();my @array_posi=();
my $read_name=''; my $site_num=0; my $real_allele;  my $heter_perce; my $posi;
my $win_snp=0;my $win_perce; my $avg_read;
my $file='';
my %event;my %read;
my $total_share=0; my $total_unique=0;
my $num_positive=0;my $num_negative=0;

$file = $depth."_event_identified_1";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
if ($b[1] eq $ID){
$event{$b[2]}="shared";
$total_share++;
}}
close I;

$file = $depth."_event_unidentified_1";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
if ($b[1] eq $ID){
$event{$b[2]}="unique";
$total_unique++;
}}
close I;

foreach my $posi (sort { $a <=> $b} keys %event) {
	if ($event{$posi} eq "shared"){
  push @win_allele,1;
	push @array_posi,$posi;
  $site_num++;
}elsif ($event{$posi} eq "unique"){
	push @win_allele,0;
	push @array_posi,$posi;
  $site_num++;
}

if (@win_allele  == $event_num){
	$num_positive=sum(@win_allele);
$win_perce=sum(@win_allele)/$event_num;

$num_negative=@win_allele-$num_positive;
my $other_posi=$total_share-$num_positive; my $other_negative=$total_unique-$num_negative;
print OUT3 "$ID\t","$site_num\t","$win_perce\t","$array_posi[0]\t","$array_posi[-1]\t","$num_positive\t","$num_negative\t","$other_posi\t","$other_negative\n";

@win_allele=(); @array_posi=();@win_read=();
}
}
$win_perce=sum(@win_allele)/@win_allele;
$num_positive=sum(@win_allele);
$num_negative=@win_allele-$num_positive;
my $other_posi=$total_share-$num_positive; my $other_negative=$total_unique-$num_negative;
print OUT3 "$ID\t","$site_num\t","$win_perce\t","$array_posi[0]\t","$array_posi[-1]\t","$num_positive\t","$num_negative\t","$other_posi\t","$other_negative\n";
}