#!/usr/bin/perl
use strict;
use warnings;

###find candidate CO reads identified with both references 
my $snp_num=10;my %tot_read;my %tot_read2;my %focal_read;my %focal_read2;my $num_read=0;my $num_read2=0; my $total_num=0;my $total_num2=0;
my $id;

my @ref=("A6","A4");
my @que=("A4","A6");
for (my $j=0; $j<@ref; $j+=1){
my $outfile1=$ref[$j]."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(OUT1, ">$outfile1");

my $file = $ref[$j]."_assembly_region_sliding_window_benchmark_pool_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		$id=$b[0];
		if (exists $tot_read{$id}){
		}else{
	$tot_read{$id}=1;
		$total_num++;
}
			if ($b[-3] == 0.5){
			if (exists $focal_read{$id}){
			}else{
		$focal_read{$id}=1;
		$num_read++;
}}
}
close FILE1;

$file = $que[$j]."_assembly_region_sliding_window_benchmark_pool_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		$id=$b[0];
		if (exists $tot_read2{$id}){
		}else{
	$tot_read2{$id}=1;
		$total_num2++;
}
			if ($b[-3] == 0.5){
				if (exists $focal_read{$id}){
			if (exists $focal_read2{$id}){
			}else{
		$focal_read2{$id}=1;
		$num_read2++;
}}}}

open(FILE,"<", "$file")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		$id=$b[0];
		if (exists $focal_read2{$id}){
			print OUT1 "$count\n";
		}}
print "$num_read\t","$num_read2\n";
}

### identify consecative parental origin switch with both references 
my $left_edge=10; my $right_edge=10; my $missing=2;
foreach my $ref (@ref){
my $outfile2=$ref."_assembly_benchmark_pool_sliding_window_breakpoint_location".$snp_num."_".$left_edge."_".$right_edge."_".$missing;
open(OUT2, ">$outfile2");

my $num_read=0; my $dup_read=0;
my %candid_read; my @loca=(); my $posi; my $id;
my $file = $ref."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);		
		$posi = $count;my $left=$b[-2]; my $right=$b[-1];
		$id = $b[0]."\t".$b[1]."\t".$b[3];
			if ($b[-3] == 0.5){
			if (exists $candid_read{$id}){
				push @loca, $posi;
				$dup_read++;
			}else{
			push @loca, $posi;
		$candid_read{$id}=1;
		$num_read++;
}}}

my $left_count=0; my $right_count=0; my $left_consist=0;my $right_consist=0;
my $read_posi; my $read_name;
for (my $i=0; $i<@loca; $i+=1){
	my @a = split("\t", $loca[$i]);
	#print "$a[0]\t","$a[1]\n";
	$left_count=0; $right_count=0; $left_consist=0;$right_consist=0;
	$read_posi=$a[5];
	$read_name = $a[0]."\t".$a[1]."\t".$a[3];

open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
$id = $b[0]."\t".$b[1]."\t".$b[3];

	if (exists $candid_read{$id}){
		my $site_posi=$b[5];
		my $prop=$b[6];
		my $name = $b[0]."\t".$b[1]."\t".$b[3];

		if ($name eq $read_name){
			if (($site_posi < $read_posi)&&($site_posi >= $read_posi-$left_edge)){
				$left_count++;
			if ($prop < 0.5){
				$left_consist=$left_consist-1;
			}elsif($prop > 0.5){
				$left_consist=$left_consist+1;
			}
}elsif (($site_posi > $read_posi)&&($site_posi <= $read_posi+$right_edge)){
		$right_count++;
	if ($prop < 0.5){
		$right_consist=$right_consist-1;
	}elsif($prop > 0.5){
		$right_consist=$right_consist+1;
	}
}}}}
close I;
print "$loca[$i]\t","$left_count\t","$right_count\t","$left_consist\t","$right_consist\n";
	if ((($left_count-abs($left_consist))<=$missing)&&(($right_count-abs($right_consist))<=$missing)){
		if ($left_consist*$right_consist<0){
		print OUT2 "$loca[$i]\t","$left_consist\t","$right_consist\n";
	}}}
}



######## identify candidate reads with double COs
my $snp_num=4;my %tot_read;my %tot_read2;my %focal_read;my %focal_read2;
my $num_snp=0; my $num_snp2=0; my $num_read=0;my $num_read2=0; my $total_num=0;my $total_num2=0;
my $inter_snp=2;
my $last_win;

for (my $j=0; $j<@ref; $j+=1){
my $file = $que[$j]."_assembly_region_sliding_window_benchmark_pool_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
		$num_snp++;
		my $region = $b[0];
		if (exists $tot_read{$region}){
		}else{
	$tot_read{$region}=1;
		$total_num++;
}
my $win_location=$b[0]."\t".$b[1]."\t".$b[2]."\t".$b[3]."\t".$b[4]."\t".$b[5];
			if ($b[-3] == 0.5){
				my @a = split ("\t", $last_win);
					if (($b[0] eq $a[0])&&($b[1] eq $a[1])&&($b[3] eq $a[3])&&(($b[5]-$a[5])>$inter_snp)){
			if (exists $focal_read{$region}){
			}else{
		$focal_read{$region}=1;
		$num_read++;
}
}
$last_win=$win_location;
}
}
close FILE1;

my $outfile3=$ref[$j]."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(OUT3, ">$outfile3");

$file = $ref[$j]."_assembly_region_sliding_window_benchmark_pool_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
    #if ($b[2] == 0.5){
		$num_snp2++;
		my $region = $b[0];
		if (exists $tot_read2{$region}){
		}else{
	$tot_read2{$region}=1;
		$total_num2++;
}
my $win_location=$b[0]."\t".$b[1]."\t".$b[2]."\t".$b[3]."\t".$b[4]."\t".$b[5];
			if ($b[-3] == 0.5){
				my @a = split ("\t", $last_win);
				if (($b[0] eq $a[0])&&($b[1] eq $a[1])&&($b[3] eq $a[3])&&(($b[5]-$a[5])>$inter_snp)){
				if (exists $focal_read{$region}){
			if (exists $focal_read2{$region}){
			}else{
		$focal_read2{$region}=1;
		$num_read2++;
}}
}
$last_win=$win_location;
}
}

open(FILE,"<", "$file")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
		my $region = $b[0];
		if (exists $focal_read2{$region}){
			print OUT3 "$count\n";
}}
}


my $left_edge=3; my $right_edge=3; my $inter_snp=2;
foreach my $ref (@ref){
my $outfile4=$ref."_assembly_gc_share_benchmark_pool_sliding_window_region_location".$snp_num."_".$left_edge."_".$right_edge;
open(OUT4, ">$outfile4");

my $num_read=0; my $dup_read=0;
my %candid_read; my @loca=(); my $posi;
my $last_win;
my $file = $ref."_share_benchmark_pool_sliding_window_prop_posi_".$snp_num.".txt";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    my @b = split("\t", $count);
my $win_location=$b[0]."\t".$b[1]."\t".$b[2]."\t".$b[3]."\t".$b[4]."\t".$b[5];
if ($b[-3] == 0.5){
	my @a = split ("\t", $last_win);
		if (($b[0] eq $a[0])&&($b[1] eq $a[1])&&($b[3] eq $a[3])&&(($b[5]-$a[5])>$inter_snp)){
		my $pair_win=$last_win.":".$win_location;
		my $id = $b[0]."\t".$b[3];
		$candid_read{$id}=1;
		push @loca, $pair_win;
$num_read++;
}
$last_win=$win_location;
}
}
print "$num_read\n";

my $left_count=0; my $right_count=0; my $left_consist=0;my $right_consist=0;my $inter_count=0;my $inter_consist=0;
for (my $i=0; $i<@loca; $i+=1){
	my @a = split(":", $loca[$i]);
	my @left_win = split("\t", $a[0]);my @right_win = split("\t", $a[1]);
	$left_count=0; $right_count=0; $left_consist=0;$right_consist=0;$inter_count=0;$inter_consist=0;

open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @b = split("\t", $count);
my $region = $b[0]."\t".$b[3];
	if (exists $candid_read{$region}){
		if ($b[0] eq $left_win[0]){
			if (($b[5] < $left_win[5])&&($b[5] >= $left_win[5]-$left_edge)){
				#if ($b[-1] < $a[-1]){
			if ($b[-3] < 0.5){
				$left_consist=$left_consist-1;
				$left_count++;
			}elsif($b[-3] > 0.5){
				$left_consist=$left_consist+1;
				$left_count++;
			}
}elsif (($b[5] > $right_win[5])&&($b[5] <= $right_win[5]+$right_edge)){
#}elsif ($b[-1] > $a[-1]){
	if ($b[-3] < 0.5){
		$right_consist=$right_consist-1;
		$right_count++;
	}elsif($b[-3] > 0.5){
		$right_consist=$right_consist+1;
		$right_count++;
	}
}elsif (($b[5] > $left_win[5])&&($b[5] < $right_win[5])){
if ($b[-3] < 0.5){
$inter_consist=$inter_consist-1;
$inter_count++;
}elsif($b[-3] > 0.5){
$inter_consist=$inter_consist+1;
$inter_count++;
}
}}}
}
close I;
if (($left_count >= $left_edge)&&($right_count >= $right_edge)&&($inter_count >= $inter_snp)){
	if ((abs($left_consist) == $left_count)&&(abs($right_consist) == $right_count)&&(abs($inter_consist) == $inter_count)){
		if (($left_consist*$right_consist>0)&&($left_consist*$inter_consist<0)){
		print OUT4 "$loca[$i]\t","$left_consist\t","$right_consist\t","$inter_consist\n";
	}}}}
}