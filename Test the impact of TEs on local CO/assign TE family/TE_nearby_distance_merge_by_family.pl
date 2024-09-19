#!/usr/bin/perl
use strict;
use warnings;

my @name=("A4","A6","A7");my @cutoff=("_10");

for (my $j=0; $j<@name; $j+=1){
foreach my $e (@cutoff){

my $i=2;my $t=0;my $c=0;
my $family=''; my $super_family=''; my $type='';
my $chr=''; my $start=0; my $end=0;


my $outfile=$name[$j]."_distance_nearbyall_type_unmerged_euchromatic_include_INE".$e;
open(OUT, ">$outfile");

my $infile2=$name[$j]."_blast_TE_family_location".$e;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);

if (($b[7]-$b[6])>=500){
	$c++;
}

if ($b[5] eq $chr){#need change $b[1] to $b[0]
my $dis=$b[6] - $end;

$start = $b[6];
$end = $b[7];
$super_family = $b[11];
$family=$b[0];
#$type= $b[5];
$i=1;

$t++;
print OUT "$dis\n";
print OUT "$count\t","$dis\t";

}else{
$t++;
if ($chr ne ''){
print OUT "NA\n";
}
print OUT "$count\t","NA\t";
$chr = $b[5];
$start = $b[6];
$end = $b[7];
$super_family = $b[11];
$family=$b[0];
#$type= $b[5];
$i=0;
}
}

print OUT "NA\n";
$t++;
print "$name[$j] ","$t ","$c\n";
}
}

my $dis=5000;
for (my $k=0; $k<@name; $k+=1){
foreach my $e (@cutoff){
my $outfile=$name[$k]."_exclu_within_".$dis."bp_distance_nearbyall_TE_euchromatic.txt";
open(OUT, ">$outfile");

my $t=0;my $c=0;
my %loca; my $family=''; my $start; my $end;

my $infile="/dfs7/grylee/yuhenh3/recombination_rate/".$name[$k]."_distance_nearbyall_type_unmerged_euchromatic_include_INE".$e;
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
$t++;
my $i=0;my $j=0;
if ($b[-2] eq "NA"){
$start=$b[6];
}elsif($b[-2] <$dis){
$i++;
}
if ($b[-1] eq "NA"){
}elsif($b[-1] <$dis){
$j++;
}

if(($i==0)&&($j==0)&&($b[0] ne "INE-1")){
print OUT "$b[0]\t","$b[5]\t","$b[6]\t","$b[7]\t","$b[11]\t","$b[-2]\t","$b[-1]\t","unmerged\n";

$c++;
}elsif(($i==0)&&($j>0)){
$start=$b[6];
$family=$b[0];
}elsif(($i>0)&&($j==0)){

if (($b[0] eq $family)&&($b[0] ne "INE-1")){
print OUT "$b[0]\t","$b[5]\t","$start\t","$b[7]\t","$b[11]\t","$b[-2]\t","$b[-1]\t","merged\n";
$c++;
}elsif ($b[0] eq "ambiguous"){
print OUT "$family\t","$b[5]\t","$start\t","$b[7]\t","$b[11]\t","$b[-2]\t","$b[-1]\t","merged\n";
$c++;
}elsif ($b[0] ne "INE-1"){
print OUT "nearby_ambiguous\t","$b[5]\t","$start\t","$b[7]\t","$b[11]\t","$b[-2]\t","$b[-1]\t","merged\n";
}
}else{
if ($b[0] eq $family){
}else{
if ($family eq "ambiguous"){
	$family = $b[0];
}elsif($b[0] eq "ambiguous"){
	$family = $family;
}else{
$family='nearby_ambiguous';
}
}}}

print "$name[$k] ","$c ","$t\n";
}}