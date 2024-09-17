#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw(sum);

my @chr=("CM010569.1","CM010570.1","CM010571.1","CM010572.1","CM010573.1");#use when ref is A6

#my @chr=("CM010576.1","CM010577.1","CM010578.1","CM010579.1","CM010580.1"); #use when ref is A7

my $ref="A6"; my $que="A4";
my $outfile1=$ref."_refer_".$que."_query_asm10_update_SNP_shared_assem_Illumia_two.vcf";

open(OUT1, ">$outfile1");
my $posi_index=0;my $chr;my $file = '';
my %ref_neu;my %que_neu;
$file=$ref."_genomic.fna";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    if($count =~ m/>/){
      my @b = split(" ", $count);
      my @c = split(">", $b[0]);
  	   $chr=$c[1];
       $posi_index=0;
}else{
    my @b = split("", $count);
		for (my $i=0; $i < @b; $i++) {
			$posi_index++;
      my $posi=$chr."\t".$posi_index;
			my $nuclotide=$b[$i];
			$ref_neu{$posi}=$nuclotide;
}}}

$file=$que."_genomic.fna";
open(FILE1,"<", "$file")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
    if($count =~ m/>/){
      my @b = split(" ", $count);
      my @c = split(">", $b[0]);
  	   $chr=$c[1];
       $posi_index=0;
}else{
    my @b = split("", $count);
		for (my $i=0; $i < @b; $i++) {
			$posi_index++;
      my $posi=$chr."\t".$posi_index;
			my $nuclotide=$b[$i];
			$que_neu{$posi}=$nuclotide;
}}}



my %ref_snp;my %que_snp;
$file = $ref."_".$ref."_genomic_sorted.g.vcf";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    if($count =~ m/#/){
}else{
    my @b = split("\t", $count);
		my $site=$b[0]."\t".$b[1];
		$ref_snp{$site}=$b[2];
	}}

	$file = $que."_".$que."_genomic_sorted.g.vcf";
	  open I, "<$file" or print "Can't open /$file\n";
		while(my $count = <I>){
	  	chomp($count);
	    if($count =~ m/#/){
	}else{
	    my @b = split("\t", $count);
			my $site=$b[0]."\t".$b[1];
			$que_snp{$site}=$b[2];
		}}


foreach my $chr (@chr){

my $total_win=0; my $window_num=1;my $SNP=0;
my $nucl=0;
$file = $ref."_refer_".$que."_query_asm10_update.var.vcf";
  open I, "<$file" or print "Can't open /$file\n";
	while(my $count = <I>){
  	chomp($count);
    if($count =~ m/#/){
		}else{
    my @b = split("\t", $count);
    if ($b[0] eq $chr){
    my @ref=split("", $b[3]);    my @alt=split("", $b[4]);
    if ((@ref == 1)&&(@alt == 1)){
      my @c=split(";", $b[7]);
      my @d= split("=", $c[0]);
      my @e= split("=", $c[1]);

		my $ref_posi = $d[1]."\t".$e[1];
    	my $que_posi = $b[0]."\t".$b[1];
		if ((exists $ref_snp{$ref_posi})||(exists $que_snp{$que_posi})){
	}else{

    if ((exists $ref_neu{$ref_posi})&&(exists $que_neu{$que_posi})){
      $nucl++;
			my $A6=uc $ref_neu{$ref_posi};
			my $A4=uc $que_neu{$que_posi};
			if (($b[3] eq $A4)&&($b[4] eq $A6)){
        $SNP++;
        print OUT1 "$count\n";
      }
  }}}}}}
#print "asm10\t","$chr\t","$nucl\t","$SNP\n";
}
