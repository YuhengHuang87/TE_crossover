#!/usr/bin/perl
use strict;
use warnings;

my %TE_neu; my $name; my $posi;
my @TE_pacbio;
my $strain="A6";
my $outfile=$strain."_TE_library_euchro.fasta";
open(OUT, ">$outfile");
my $outfile1=$strain."_TE_euchro.out";
open(OUT1, ">$outfile1");
my $file_name=$strain."_genomic.fna";
my $infile1=$file_name.".out";
open(FILE1,"<", "$infile1")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		my @b = split(" ", $count);
		my $chr = $b[4]; my $left = $b[5]; my $right = $b[6];

		#if ((($chr eq "CM010541.1")&&($left>206715)&&($right<18899017))||(($chr eq "CM010542.1")&&($left>123450)&&($right<19558299))||(($chr eq "CM010543.1")&&($left>8597017)&&($right<24417585))||(($chr eq "CM010544.1")&&($left>173141)&&($right<18439213))||(($chr eq "CM010545.1")&&($left>10069791)&&($right<32445046))){ #Dmel for A4
		if ((($chr eq "CM010569.1")&&($left>144265)&&($right<18689904))||(($chr eq "CM010570.1")&&($left>251300)&&($right<19760175))||(($chr eq "CM010571.1")&&($left>8030989)&&($right<23829870))||(($chr eq "CM010572.1")&&($left>164091)&&($right<18449201))||(($chr eq "CM010573.1")&&($left>8729782)&&($right<31194759))){ #Dmel for A6
		#if ((($chr eq "CM010576.1")&&($left>204176)&&($right<18599896))||(($chr eq "CM010577.1")&&($left>134798)&&($right<19487247))||(($chr eq "CM010578.1")&&($left>9035199)&&($right<24726584))||(($chr eq "CM010579.1")&&($left>292296)&&($right<18448828))||(($chr eq "CM010580.1")&&($left>9655508)&&($right<32200062))){ #Dmel for A7

		if (($b[10] =~ /LINE/)||($b[10] =~ /LTR/)||($b[10] =~ /DNA/)||($b[10] =~ /Unknown/)){
	 $name=$b[10].":".$chr.":".$left.":".$right;
push @TE_pacbio,$name;
print OUT1 "$count\n";
	for (my $i=$left; $i<=$right; $i+=1){
		$posi = $chr."\t".$i;
		$TE_neu{$posi}=$name;
}}
}
}

my $num_fam=scalar(@TE_pacbio);
print "$num_fam\n";

my %TE_seq; my $chr;my $posi_index=0;
my $infile2=$file_name;
open(FILE1,"<", "$infile2")||die"$!";
	while(my $count = <FILE1>){
		chomp($count);
		if($count =~ m/>/){
			my @a = split(" ", $count);
			my @c = split(">", $a[0]);
			$chr = $c[1];
			 $posi_index=0;
	}else{
		my @d = split("", $count);
		for (my $i=0; $i < @d; $i++) {
			$posi_index++;
			my $site=$chr."\t".$posi_index;
			if (exists $TE_neu{$site}){
				if (exists $TE_seq{$TE_neu{$site}}){
					$TE_seq{$TE_neu{$site}}=$TE_seq{$TE_neu{$site}}.$d[$i];
				}else{
					$TE_seq{$TE_neu{$site}}=$d[$i];
				}
			}
	}}}


foreach my $key (keys %TE_seq) {
	print OUT ">","$key\n","$TE_seq{$key}\n";
}
