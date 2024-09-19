#!/usr/bin/perl
use strict;
use warnings;

my @chromo=("CM010569.1","CM010570.1","CM010571.1","CM010572.1","CM010573.1");
my @A4_chr=("CM010541.1","CM010542.1","CM010543.1","CM010544.1","CM010545.1");
my @chr=("X","2L","2R","3L","3R");


for (my $i=0; $i < scalar(@chromo); $i++) {
my $chromo=$chromo[$i];
my $A4_chr=$A4_chr[$i];

my $outfile1=$chromo."_breakpoint_read_que_TE_within_prop_posi";
open(OUT1, ">$outfile1");

my $outfile2=$chromo."_breakpoint_read_que_TE_notwithin_prop_posi";
open(OUT2, ">$outfile2");

my %map_TE;
my $infile1="A6_exclu_within_5000bp_distance_nearbyall_TE_euchromatic.txt"; #TE profile generated from scripts in folder /Test the impact of TEs on local CO/
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
		my $TE_length=$b[3]-$b[2]+1;
		if ($b[1] eq $chromo){
            my $loca = $b[1]."\t".$b[2]."\t".$b[3];
            $map_TE{$loca}=1;
        }}

my $file;
my %posi_read;
my $start; my $end; my $frag_id; my $id; my $n=0;
$file= "A6_assembly_0.1_sliding_window_".$chromo."_prop_posi_10.txt";#generated from the CO calling pipeline 
open(FILE,"<", "$file")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
        if ($b[1] eq $chromo){
        $id = $b[0]."\t".$b[2];
        if ($n == 0){
            $frag_id=$id; $start=$b[-2];
        }else{
        if ($id eq $frag_id){
            $end = $b[-1];
        }else{
		$posi_read{$frag_id}=$chromo."\t".$start."\t".$end;
        $frag_id = $id; $start=$b[-2];
        }}
        $n++;
        }}
        $posi_read{$frag_id}=$chromo."\t",$start."\t".$end;
#print "$n\n";
close FILE;

my $hit=0;
foreach my $match (keys %posi_read){
    print "$posi_read{$match}\n";
    $hit=0;
		my @range = split("\t", $posi_read{$match});
        foreach my $loca (keys %map_TE){
            my @b = split("\t", $loca);
		if ($b[0] eq $range[0]){
                if (($b[2]<$range[1])||($b[1]>$range[2])){
                }else{
                    $hit++;
                }}}
  if ($hit > 0){
    print OUT1 "$match\t","$posi_read{$match}\n";
  }else{
    print OUT2 "$match\t","$posi_read{$match}\n";
  }}
}


for (my $i=0; $i < scalar(@chromo); $i++) {
my $chr=$chr[$i];
my $chromo=$chromo[$i];
my $A4_chr=$A4_chr[$i];

my $outfile1=$chromo."_breakpoint_read_que_SV_within_prop_posi";
open(OUT1, ">$outfile1");

my $outfile2=$chromo."_breakpoint_read_que_SV_notwithin_prop_posi";
open(OUT2, ">$outfile2");

my %map_TE; 
my $infile1="master_table_d10_100bp.txt";#obtained from Chakraborty et al. 2019
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split("\t", $count);
            if ($b[7] =~ m/ME/){}else{
            my $left=0; my $right=0;
        if($b[8] =~ m/\;/){
            #if(($b[8] =~ m/a6/)&&($b[8] =~ m/a7/)){$n++;}else{
            my @d = split(/\;/,$b[8]);
            my @left = split(/\;/,$b[9]);my @right = split(/\;/,$b[10]);
            #print "@left\n";
            for (my $j=0; $j < @d; $j++) {
        my @a = split(/\./,$d[$j]);
        #print "$a[0]\n";
        if ($a[0] eq "a6"){
		if ($a[1] eq $chr){
            my $loca = $chromo."\t".$left[$j]."\t".$right[$j];
            print "$loca\n";
            $map_TE{$loca}=1;
        }}}
        #}
        }else{
        my @a = split(/\./,$b[8]);
        if ($a[0] eq "a6"){
		if ($a[1] eq $chr){
            my $loca = $chromo."\t".$b[9]."\t".$b[10];
            $map_TE{$loca}=1;
            #print "$loca\n";
        }
		}
        }}}

my $file;
my %posi_read;
my $start; my $end; my $frag_id; my $id; my $n=0;
$file= "A6_assembly_0.1_sliding_window_".$chromo."_prop_posi_10.txt";#generated from the CO calling pipeline
open(FILE,"<", "$file")||die"$!";
	while(my $count = <FILE>){
		chomp($count);
		my @b = split("\t", $count);
        if ($b[1] eq $chromo){
        $id = $b[0]."\t".$b[2];
        if ($n == 0){
            $frag_id=$id; $start=$b[-2];
        }else{
        if ($id eq $frag_id){
            $end = $b[-1];
        }else{
		$posi_read{$frag_id}=$chromo."\t".$start."\t".$end;
        $frag_id = $id; $start=$b[-2];
        }}
        $n++;
        }}
        $posi_read{$frag_id}=$chromo."\t",$start."\t".$end;
#print "$n\n";
close FILE;

my $hit=0;
foreach my $match (keys %posi_read){
    $hit=0;
		my @range = split("\t", $posi_read{$match});
        foreach my $loca (keys %map_TE){
            my @b = split("\t", $loca);
		if ($b[0] eq $range[0]){
                if (($b[2]<$range[1])||($b[1]>$range[2])){
                }else{
                    $hit++;
                }}}
  if ($hit > 0){
    print OUT1 "$match\t","$posi_read{$match}\n";
  }else{
    print OUT2 "$match\t","$posi_read{$match}\n";
  }}
}