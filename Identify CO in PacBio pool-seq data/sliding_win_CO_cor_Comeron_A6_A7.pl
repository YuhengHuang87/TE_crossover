#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);

my $snp_num=10;
my $outfile="cor_Comeron_A6_A7_rate_depth_normalized.txt";
open(OUT, ">$outfile");


my $read_name=''; my $site_num=0; my $real_allele; my $win_perce; my $heter_perce; 
my $posi;my $homo_rec; my $last_region;
my $win_snp=0; my $file='';

my @chr=("2L","2R","3L","3R","X");
foreach my $chr (@chr){
my %array_posi;
$file = "comeron_pool_A6_A7_100kb_rate_depth_normalized";
open I, "<$file" or print "Can't open /$file\n";
while(my $count = <I>){
chomp($count);
my @a = split(' ', $count);

if ($a[0] eq $chr){
  $posi = $a[1];
  $homo_rec=$a[2]."\t".$a[6]."\t".$a[11];
    $array_posi{$posi}=$homo_rec;
}
}

my @comeron; my @A6; my @A7;
foreach my $key (sort { $a <=> $b} keys %array_posi) {
    my @b = split("\t", $array_posi{$key});
    push @comeron, $b[0];
    push @A6, $b[1];
    push @A7, $b[2];

if (@comeron  == $snp_num){
my $win_comeron=sum(@comeron);
my $win_A6=sum(@A6);
my $win_A7=sum(@A7);

print OUT "$chr\t","$key\t","$win_comeron\t","$win_A6\t","$win_A7\n";
splice @comeron,0,1;
splice @A6,0,1;
splice @A7,0,1;
}
}
close I;
}