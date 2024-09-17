#!/usr/bin/perl
use strict;
use warnings;
my $snp_num=10;my $left_edge=10; my $right_edge=10;  my $missing=2; 
my $ref="A6";

#### filter out reads show more than one CO
my $outfile=$ref."_benchmark_pool_CO_location_both_parent_exclu_gc_".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open(OUT, ">$outfile");

my $lin=0; my $posi_A6_1;my $posi_A6_2;my $lin_2=0;my $lin_A4=0; my $lin_A6=0; my $co_gc_share=0;
my %A4_gc_read;my %gc;my $id;

my $infile="A4_assembly_gc_share_benchmark_pool_sliding_window_region_location4_3_3";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		$id = $b[0];
		$A4_gc_read{$id}=1;
}
$infile="A6_assembly_gc_share_benchmark_pool_sliding_window_region_location4_3_3";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
	if (exists $A4_gc_read{$b[0]}){
		$id = $b[0]."\t".$b[3];
		$gc{$id}=1;
		$lin++;
}}

my %A4_co_read; my %unique_A6_reads;
$infile="A4_assembly_benchmark_pool_sliding_window_breakpoint_location10_10_10_2";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		$id = $b[0];
		if (exists $A4_co_read{$id}){
		}else{
		$A4_co_read{$id}=1;
		$lin_A4++;
}}

$infile="A6_assembly_benchmark_pool_sliding_window_breakpoint_location10_10_10_2";
open(FILE,"<", "$infile")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		$id = $b[0]."\t".$b[3];	
	if (exists $A4_co_read{$b[0]}){
			print OUT "$count\n";
			if (exists $unique_A6_reads{$id}){}else{
						$lin_2++;
						$unique_A6_reads{$id}=1;
			if (exists $gc{$id}){
		$co_gc_share++;
				}else{
			}}}}
print "$lin\t","$lin_2\t","$lin_A4\t","$co_gc_share\n";

#### filter out reads where the CO location was within 2kb from the edge of alignment
my $dis_cutoff=2000;#can change the distance of CO to the edge of alignment 
my $num_read=0;my $num_read2=0;
my $outfile1=$ref."_CO_shared_each_read_fragment_distance_edge_".$dis_cutoff."_".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open(OUT1, ">$outfile1");

my %posi_read;my %focal_read; my $id;
my %A4_co_read;
my $file="A4_assembly_benchmark_pool_sliding_window_breakpoint_location".$snp_num."_10_10_".$missing;
open(FILE,"<", "$file")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		$id = $b[0]."\t".$b[3];
		$A4_co_read{$id}=($b[-3]+$b[-4])/2;;
}

my %A4_posi;
foreach my $read_name (keys %A4_co_read){
	#print "$read_name\n";
	my %focal_read3;
	my @a = split("\t", $read_name);
	my $left=0;my $right=0;
	$file = "A4_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
		while(my $count = <FILE1>){
			chomp($count);
	    my @b = split("\t", $count);
			$id = $b[0]."\t".$b[3];
	    if ($id eq $read_name){
				if (exists $focal_read3{$id}){
				}else{
			$focal_read3{$id}=1;
			$left=($b[-2]+$b[-1])/2;
	}
	$right=($b[-2]+$b[-1])/2;
	}
	}
	my $dis1=$A4_co_read{$read_name}-$left; my $dis2=$right-$A4_co_read{$read_name};
	if (($dis1>=$dis_cutoff)&&($dis2>=$dis_cutoff)){
		$A4_posi{$a[0]}=1;
	}
}

$file=$ref."_benchmark_pool_CO_location_both_parent_exclu_gc_".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
	my $chr=$b[1]; my $left=$b[-4]; my $right=$b[-3]; 
	if ((($chr eq "CM010569.1")&&($left>144265)&&($right<18689904))||(($chr eq "CM010570.1")&&($left>251300)&&($right<19760175))||(($chr eq "CM010571.1")&&($left>8030989)&&($right<23829870))||(($chr eq "CM010572.1")&&($left>164091)&&($right<18449201))||(($chr eq "CM010573.1")&&($left>8729782)&&($right<31194759))){ #Dmel for A6
		$id = $b[0]."\t".$b[3];
		if (exists $focal_read{$id}){}else{
		$posi_read{$id}=($b[-3]+$b[-4])/2;
		$focal_read{$id}=$count;
		$num_read++;
}}}
close FILE1;

foreach my $read_name (keys %focal_read){
my $left=0;my $right=0;my %focal_read2;
my @a = split("\t", $read_name);

$file = $ref."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		$id = $b[0]."\t".$b[3];
		if ($id eq $read_name){
		if (exists $focal_read2{$id}){
			}else{
		$focal_read2{$id}=1;
		$left=($b[-2]+$b[-1])/2;
}
$right=($b[-2]+$b[-1])/2;
}}
my $dis1=$posi_read{$read_name}-$left; my $dis2=$right-$posi_read{$read_name};
if (($dis1>=$dis_cutoff)&&($dis2>=$dis_cutoff)){
	if (exists $A4_posi{$a[0]}){
	print OUT1 "$focal_read{$read_name}\n";
	$num_read2++;
}}}


#### further zoom in the location of CO
my $outfile2=$ref."_CO_share_singlePool_sliding_window_location_zoom_in_".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open(OUT2, ">$outfile2");

my $file='';
my @breakpoint=();
my $indiv;
$file=$ref."_CO_shared_each_read_fragment_distance_edge_".$dis_cutoff."_".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
my @b = split("\t", $count);
$tot++;
	$crossover++;
	my $range=$b[0]."\t".$b[1]."\t".$b[-4]."\t".$b[-3];
	push @breakpoint, $range;
}
	my $front_one='';	 my $front_two=''; my $win_num=0;
	my $left_bound='';my $right_bound='';
for (my $i=0; $i<@breakpoint; $i+=1){
	my @a = split("\t", $breakpoint[$i]);
	$left_bound='';$right_bound='';

$file = $ref."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
my @b = split("\t", $count);
if ($b[0] eq $a[0]){
	if ($b[1] eq $a[1]){
		if(($b[2] >= $a[2])&&($b[2] <= $a[3])){
			$win_num++;
		}
	if ($win_num == 5){
		$left_bound=$b[2];
}elsif($win_num == 6){
		$right_bound=$b[2];
	}
}
}
}
print OUT2 "$a[0]\t","$a[1]\t","$left_bound\t","$right_bound\n";
	$win_num=0;
}
print "$tot\t","$crossover\t","$final\n";