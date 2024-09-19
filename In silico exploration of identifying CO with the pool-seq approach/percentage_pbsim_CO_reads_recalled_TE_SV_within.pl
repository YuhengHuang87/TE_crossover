#!/usr/bin/perl
use strict;
use warnings;

my $m=0; my $n=0;my $p=0; my $q=0;

#my $outfile=$r."_breakpoint_read_identified_contain_TE_orNot.txt";
my $outfile=$r."_breakpoint_read_identified_contain_SV_orNot.txt";
open(OUT, ">$outfile");

my %identified; my %unided; 
my $infile1="A6_co_gc_shared_each_read_fragment_cutoff_2000_10_10_10_2";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
my $id = $b[0]."\t".$b[1];
	$identified{$id}=$count;
		}
close FILE1;

my @chromo=("CM010569.1","CM010570.1","CM010571.1","CM010572.1","CM010573.1");
foreach my $chromo (@chromo){

#my $infile2=$chromo."_breakpoint_read_que_TE_within_prop_posi";
my $infile2=$chromo."_breakpoint_read_que_SV_within_prop_posi";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		my @b = split("\t", $count);
		my $length = $b[-1]-$b[-2]+1;
        my $id = $b[0]."\t".$b[2];
if (exists $identified{$id}){
	$m++;
	print OUT "ID\t","conTE\t","$length\n";
}else{
	$n++;
	print OUT "unID\t","conTE\t","$length\n";
}
}

#$infile2=$chromo."_breakpoint_read_que_TE_notwithin_prop_posi";
$infile2=$chromo."_breakpoint_read_que_SV_notwithin_prop_posi";

open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		my @b = split("\t", $count);
my $length = $b[-1]-$b[-2]+1;
        my $id = $b[0]."\t".$b[2];
        if (exists $identified{$id}){
	$p++;
	print OUT "ID\t","notTE\t","$length\n";
}else{
	$q++;
	print OUT "unID\t","notTE\t","$length\n";
}
}
print "$m\t","$n\t","$p\t","$q\n";
}