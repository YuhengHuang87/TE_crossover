#!/usr/bin/perl
use strict;
use warnings;

my @name=("A4","A6","A7");my @cutoff=("10");
foreach my $name (@name){
foreach my $e (@cutoff){
my $t=0;my $c=0;
my %fam_id;
my $outfile=$name."_blast_TE_euchro_canonical_unique".$e;
open(OUT, ">$outfile");
my $outfile1=$name."_blast_TE_euchro_canonical_multi".$e;
open(OUT1, ">$outfile1");

my $infile1="dmel-all-transposon-r6.32_canonical";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	if ($count =~ />/){
		my @b = split(" ", $count);
	my @a = split(">", $b[0]);
	if ($a[1] =~ /;/){
	my @c = split(";", $a[1]);
	my $fb_ID=$c[0];
	$fam_id{$fb_ID}=$fb_ID;
$c++;
}

}}

$infile1=$name."_TE_library_euchro.fasta";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	if ($count =~ />/){
		$t++;
}}


my %te; my %fam_count;
my $infile2=$name."_combineTE_canonical_db_80_blast_euchro_".$e;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	my @b = split(" ", $count);
	#print "$b[1]\n";
	my $super_family=$b[0];
if ($b[1] =~ /;/){
my @d=split(";", $b[1]);
my $name=$d[0];
if (exists $te{$super_family}){
if ($fam_count{$super_family}== 1){
if ($te{$super_family} ne $fam_id{$name}){
$te{$super_family} = $te{$super_family}." ".$fam_id{$name};
$fam_count{$super_family}=$fam_count{$super_family}+1;
}
}else{
my @ele=split(" ", $te{$super_family});
my $k=0;
foreach my $ele (@ele){
if ($ele eq $fam_id{$name}){
$k++;
}}
if ($k == 0){
$te{$super_family} = $te{$super_family}." ".$fam_id{$name};
$fam_count{$super_family}=$fam_count{$super_family}+1;
}
}
}else{
$te{$super_family} = $fam_id{$name};
$fam_count{$super_family}= 1;
}
#print OUT "$fam_id{$b[1]} ","$count\n";
}}

my $i=0;my $j=0;
foreach my $key (sort keys %te){
$i++;
if ($fam_count{$key}==1){
print OUT "$key ","$fam_count{$key}\t","$te{$key}\n";
}else{
print OUT1 "$key ","$fam_count{$key}\t","$te{$key}\n";
$j++;
}}
print "$e ","$name ","$t ","$i ","$j\n";
}
}

foreach my $name (@name){
foreach my $e (@cutoff){
my $t=0;

my $outfile=$name."_blast_TE_family_unique_euchro_".$e;
open(OUT, ">$outfile");
my $outfile1=$name."_blast_TE_family_multi_euchro_".$e;
open(OUT1, ">$outfile1");

my %fam_id;
my $infile1="dmel-all-transposon-r6.32";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	if ($count =~ />/){
		my @b = split(" ", $count);
	my @a = split(">", $b[0]);
	my $fb_ID=$a[1];
	my @c = split("=", $b[3]);
	my @d = split(/\{/, $c[1]);
	my $family_name=$d[0];
	$fam_id{$fb_ID}=$family_name;
$t++;
}}
my %can_id; my $uniq_can=0;

$infile1=$name."_blast_TE_euchro_canonical_unique".$e;
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
	$can_id{$b[0]}=$b[2];
$uniq_can++;
}

my %te; my %fam_count;
my $infile2=$name."_combineTE_canonical_db_80_blast_euchro_".$e;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	my @b = split(" ", $count);
	#print "$b[1]\n";
	if ($b[1] =~ /;/){
}else{
	my $super_family=$b[0];
if (exists $can_id{$super_family}){
}else{
my $map_family;
if ($b[1] =~/FBti/){
$map_family=$fam_id{$b[1]};
}else{
my @d = split("_",$b[1]);
$map_family=$d[-1];
}
if ($map_family ne "ambiguous"){
if (exists $te{$super_family}){
if ($fam_count{$super_family}== 1){

if ($te{$super_family} ne $map_family){
$te{$super_family} = $te{$super_family}." ".$map_family;
$fam_count{$super_family}=$fam_count{$super_family}+1;
}
}else{
my @ele=split(" ", $te{$super_family});
my $k=0;
foreach my $ele (@ele){
if ($ele eq $map_family){
$k++;
}}
if ($k == 0){
$te{$super_family} = $te{$super_family}." ".$map_family;
$fam_count{$super_family}=$fam_count{$super_family}+1;
}
}
}else{
$te{$super_family} = $map_family;
$fam_count{$super_family}= 1;
}}
#print OUT "$fam_id{$b[1]} ","$count\n";
}
}}

my $i=0;my $j=0;
foreach my $key (sort keys %te){
if ($fam_count{$key}==1){
print OUT "$key ","$fam_count{$key}\t","$te{$key}\n";
$i++;
}else{
print OUT1 "$key ","$fam_count{$key}\t","$te{$key}\n";
$j++;
}}

print "$name ","$uniq_can ","$i ","$j\n";
}
}

for (my $j=0; $j<@name; $j+=1){
foreach my $e (@cutoff){
my $t=0;
my $outfile=$name[$j]."_blast_TE_family_multi_can_select_based_on_length_".$e;
open(OUT, ">$outfile");
my %fam_id;
my $infile1="dmel-all-transposon-r6.32";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	if ($count =~ />/){
		my @b = split(" ", $count);
#print "$b[0]\n";
	my @a = split(">", $b[0]);
	my $fb_ID=$a[1];
	my @c = split("=", $b[3]);
	my @d = split(/\{/, $c[1]);
	my $family_name=$d[0];
	$fam_id{$fb_ID}=$family_name;
$t++;
}
}

my %multi_id; my $multi_num=0;
$infile1=$name[$j]."_blast_TE_family_multi_euchro_".$e;
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
	$multi_id{$b[0]}=$count;
$multi_num++;
}

my $identified=0; my $ambiguous=0;
foreach my $key (sort keys %multi_id){

my %leng; my %fam_count;
my $infile2=$name[$j]."_combineTE_canonical_db_80_blast_euchro_".$e;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
	my @b = split(" ", $count);
	if ($b[1] =~ /;/){
}else{

	my $super_family=$b[0];
if ($super_family eq $key){

my $map_family;
if ($b[1] =~/FBti/){
$map_family=$fam_id{$b[1]};
}else{
my @d = split("_",$b[1]);
$map_family=$d[-1];
}

if (exists $leng{$map_family}){
	if ($b[3] > $leng{$map_family}){
$leng{$map_family}=$b[3];
	}else{
		$leng{$map_family}=$leng{$map_family};
	}
}else{
$leng{$map_family}=$b[3];
}
}}
}

my $total_len=0; my @id=();
foreach my $key2 (sort { $leng{$b} <=> $leng{$a} } keys %leng) {
$total_len+=$leng{$key2};
push @id, $key2;
}

my $prop1=$leng{$id[0]}/$total_len;
my $prop2=$leng{$id[1]}/$total_len;

if ($prop1>0.5){
print OUT "$key ","$id[0] ","$prop1\n";
$identified++;
}else{
print OUT "$key ","ambiguous ","$id[0] ","$id[1] ","$prop1 ","$prop2\n";
$ambiguous++;
}
}
print "$name[$j] ","$identified ","$ambiguous\n";
}
}

for (my $j=0; $j<@name; $j+=1){
foreach my $e (@cutoff){

my $outfile=$name[$j]."_blast_TE_family_location_".$e;
open(OUT, ">$outfile");

my $outfile2=$name[$j]."_blast_TE_family_location_prec_map_".$e;
open(OUT2, ">$outfile2");

my %te_id; my $c=0;my $t=0;
my $infile;

$infile=$name[$j]."_blast_TE_euchro_canonical_unique".$e;
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
	$te_id{$b[0]}=$b[2];
$c++;
}

$infile=$name[$j]."_blast_TE_family_unique_euchro_".$e;
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
	$te_id{$b[0]}=$b[2];
$c++;
}

$infile=$name[$j]."_blast_TE_family_multi_can_select_based_on_length_".$e;
open(FILE1,"<", "$infile")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
	$te_id{$b[0]}=$b[1];
$c++;
}

my %super_num; my %super_map;
my $infile1=$name[$j]."_TE_euchro.out";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		my $range= $b[10].":".$b[4].":".$b[5].":".$b[6];

$t++;
my $id=$b[10];
if (exists $super_num{$id}){
$super_num{$id}=$super_num{$id}+1;
}else{
$super_num{$id}=1;
}

if (exists $te_id{$range}){
print OUT "$te_id{$range} ","$count\n";
if (exists $super_map{$id}){
$super_map{$id}=$super_map{$id}+1;
}else{
$super_map{$id}=1;
}}}
my %perc;
foreach my $key (keys %super_num){
my $r;
if (exists $super_map{$key}){
$r=$super_map{$key}/$super_num{$key};
}else{
$r=0;
}
$perc{$key}=$r;
}

foreach my $key2 (sort { $perc{$a} <=> $perc{$b} } keys %perc) {
print OUT2 "$key2\t","$super_num{$key2}\t","$perc{$key2}\n";
}


print "$name[$j]\t","$c\t","$t\n";
}
}