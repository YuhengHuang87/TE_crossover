#!/usr/bin/perl
use strict;
use warnings;

my $chromo="CM010569.1";
my $depth=15;
my $outfile=$chromo."_CO_read_".$depth."x.fastq";
open(OUT, ">$outfile");

my $outfile1=$chromo."_CO_read_".$depth."x_location";
open(OUT1, ">$outfile1");

my @co_location_individual_chrom;
my $infile="A6_A4_recombine_individual_arms_".$chromo."_location";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		my $is_odd = $b[2] % 2;
		if ($is_odd==0){
push @co_location_individual_chrom, $b[2];
}else{
	push @co_location_individual_chrom, $b[4];
}
}

for (my $i=1; $i<=400; $i++){
	my %read; my $index;
	my @a = split("", $i);
	if (@a==1){
		$index="000".$i;
	}elsif(@a==2){
		$index="00".$i;
	}elsif(@a==3){
		$index="0".$i;
	}elsif(@a==4){
		$index=$i;
	}
my $j=0; my $loca=$co_location_individual_chrom[$i-1];
my $read_loca;
my $infile1=$chromo."_".$depth."_".$index.".maf";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		if($count eq "a"){
		$j=0;
	}elsif($count =~ m/ref/){
		if (($loca > $b[2]+$edge)&&($loca < $b[2]+$b[3]-$edge)){
			$j=1;
			$read_loca=$loca."\t".$b[2]."\t".$b[3];
		}
	}elsif($count =~ m/s/){
		if ($j==1){
			$read{$b[1]}=1;
			print OUT1 "$chromo\t","$b[1]\t","$read_loca\n";
		}
	}
}
close FILE1;

my $k=0;
my $infile2=$chromo."_".$depth."_".$index.".fastq";
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
  	chomp($count);
		if($count =~ m/@/){
my @b = split("@", $count);
if (exists $read{$b[1]}){
	$k=1;
print OUT "$count","_","$chromo\n";
}else{
	$k=0;
}
}elsif($count =~ m/S/){
	if ($k==1){
		print OUT "$count","_","$chromo\n";
	}
}else{
	if ($k==1){
		print OUT "$count\n";
	}
}
}
}
#+S1_4579_CM010569.1