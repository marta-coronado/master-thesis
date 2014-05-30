#!/usr/bin/perl

use Getopt::Long;

GetOptions(
  'coding=s' => \$file1,           # GENE
  'shortintrons=s' => \$file2,     # GENE.SHORTINTRON.65_8-30
  'filter=s' => \$file3,           # G2 (pero lee una línea de los filtros que le pasamos con bash) 
);

open(INPUT1, "$file1") or die "Error 1";
@archiu1 = <INPUT1>;	
shift(@archiu1);

open(INPUT2, "$file2") or die "Error 2";
@archiu2 = <INPUT2>;
shift(@archiu2);

# open(INPUT3, "$file3") or die "Error 3";
# @archiu3 = <INPUT3>;

for ($cnt = 63; $cnt >= 0; $cnt--) # to sum the SFS
{
	unshift(@col,$cnt);
	@col = sort {$a <=> $b} @col;
}
	
my @sum0 = ("3","4","5","6"); #for divergence.txt 0fold                    mdmel_0f	  mdyak_0f	  seg_0f   div_0f	
my @sum2 = ("69","70","71","72"); #for divergence.txt 2fold                mdmel_2f   mdyak_2f    seg_2f   div_2f
my @sum_neutral = ("93","94","95","96"); #for divergence.txt short introns mdmel_ins  mdyak_ins   seg_ins  div_ins
# my @sum_neutral = ("46","47","48","49"); #for divergence.txt 4fold       
# print "@sum_neutral\n";

my @sfs0 = ("25");  #for sfs.txt 0fold                 sfs_0f
my @sfs2 = ("91"); #for sfs.txt 2fold                  sfs_2f
my @sfs_neutral = ("97"); #for sfs.txt short introns   sfs_ins
# my @sfs_neutral = ("68"); #for sfs.txt 4fold 
# print "@sfs_neutral\n";


### FILTERS ###

chomp($file3); # eliminar intro del final de línea

@FILE1 = split(/:/,$file3); # separar campos por ':'


$chr = $FILE1[0];            # chr
$mutation1 = $FILE1[1];      # mutation
$mutation2 = $FILE1[2];      # mutation
$recombination1 = $FILE1[3]; # recombination
$recombination2 = $FILE1[4]; # recombination 
$domain = $FILE1[5];
$m_complexity1 = $FILE1[6];
$m_complexity2 = $FILE1[7];

#$watterson1 = $FILE1[5];     # watterson
#$watterson2 = $FILE1[6];     # watterson

=xx
$var = 1;
$cal = 0;
$dom = 0;
$lines = 128;
while ($var < $lines)
{
	$cal = 1/$var;
	$dom = $dom+$cal;
	$var++;
}
=cut
#$chr = $FILE1[0];          # chr
#$breadth1 = $FILE1[1];     # bias_exp
#$breadth2 = $FILE1[2];     # bias_exp
#$expression1 = $FILE1[3];  # max_exp
#$expression2 = $FILE1[4];  # max_exp
#$domain = $FILE1[5];       # chromatin domain

#print "\n\n  chr:\t\t".$chr."\n";
#print "  breadth:\t".$breadth1."-".$breadth2."\n";
#print "  expression:\t".$expression1."-".$expression2."\n";
#print "  domain:\t".$domain."\n\n\n";


=xx
$recombination1 = $FILE1[1];
$recombination2 = $FILE1[2]; 
$mutation1 = $FILE1[3];  
$mutation2 = $FILE1[4];
$state = $FILE1[5];
$transcripts1 = $FILE1[5];
$transcripts2 = $FILE1[6];
$size1 = $FILE1[7];
$size2 = $FILE1[8];
$distance1 = $FILE1[9];
$distance2 = $FILE1[10];
$exons1 = $FILE1[11];
$exons2 = $FILE1[12];
$breadth1 = $FILE1[13];
$breadth2 = $FILE1[14];
$expression1 = $FILE1[15];
$expression2 = $FILE1[16];
$m_complexity1 = $FILE1[17];
$m_complexity2 = $FILE1[18];
print "$chr\t$mutation1\n";




my $m = $FILE2[3] + $FILE2[46] + $FILE2[69];
and $m > $size1 and $m <= $size2
and $FILE2[31] > $breadth1 and $FILE2[31] <= $breadth2
and $FILE2[28] > $transcripts1 and $FILE2[28] <= $transcripts2
and $FILE2[30] > $expression1 and $FILE2[30] <= $expression2 
my $complex = $FILE2[28]/$FILE2[27] if ($FILE2[27] > 0);
and $complex > $m_complexity1 and $complex <= $m_complexity2
and $FILE2[100] eq $domain 
and $FILE2[29] > $distance1 and $FILE2[29] <= $distance2 
and $FILE2[27] > $exons1 and $FILE2[27] <= $exons2 
  
=cut

### 0-FOLD ###

foreach $col_sum(@sum0) 
{
	my $results = 0;
	foreach $B(@archiu1) 
	{
		chomp($B);
		my @FILE2 = split(/\t/,$B);
		my $div = $FILE2[96]/$FILE2[94]; # ins
		my $complex = $FILE2[28]/$FILE2[27] if ($FILE2[27] > 0);
		#my $m = $FILE2[3] + $FILE2[46] + $FILE2[69];
		#my $K_wat = $FILE2[95]/$FILE2[93]; 
		#my $watterson = $K_wat/$dom;
		if ($chr eq $FILE2[0] and $FILE2[41] > $recombination1 and $FILE2[41] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $FILE2[100] eq $domain and $complex > $m_complexity1 and $complex <= $m_complexity2) 
		{
			$results += $FILE2[$col_sum]; # breadth, expression, chromatin domain
		}
	}
	print $results." ";
}

foreach $col_sfs (@sfs0) 
{
	my @site_frq = "";
	my @site_frq_cell = "";
	foreach $cell(@col) 
	{
		my $results = 0;
		foreach $B(@archiu1) 
		{
			chomp($B);
			my @FILE2 = split(/\t/,$B);
			my $div = $FILE2[96]/$FILE2[94]; # ins
			my $complex = $FILE2[28]/$FILE2[27] if ($FILE2[27] > 0);
			#my $m = $FILE2[3] + $FILE2[46] + $FILE2[69];
			#my $K_wat = $FILE2[95]/$FILE2[93]; 
			#my $watterson = $K_wat/$dom;
			if ($chr eq $FILE2[0] and $FILE2[41] > $recombination1 and $FILE2[41] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $FILE2[100] eq $domain and $complex > $m_complexity1 and $complex <= $m_complexity2)		
			{
				my @site_frq = $FILE2[$col_sfs];
				my @site_frq_cell = split(/:/,$site_frq[0]);
				$results += $site_frq_cell[$cell];
			}
		}
		print $results.":";
	}
	print "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0 ";
}


### SHORT INTRONS or 4 fold sites ###

foreach $col_sum(@sum_neutral) 
{
	my $results = 0;
	foreach $B(@archiu2) 
	{
		chomp $B;
		my @FILE2 = split(/\t/,$B);
		my $div = $FILE2[96]/$FILE2[94]; # ins
		my $complex = $FILE2[28]/$FILE2[27] if ($FILE2[27] > 0);
		#my $m = $FILE2[3] + $FILE2[46] + $FILE2[69];
		#my $K_wat = $FILE2[95]/$FILE2[93]; 
		#my $watterson = $K_wat/$dom;
		if ($chr eq $FILE2[0] and $FILE2[41] > $recombination1 and $FILE2[41] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $FILE2[100] eq $domain and $complex > $m_complexity1 and $complex <= $m_complexity2)
		{
			$results += $FILE2[$col_sum] ;
		} 
	}
	print $results." ";
}

foreach $col_sfs(@sfs_neutral) 
{
	my @site_frq = "";
	my @site_frq_cell = "";
	foreach $cell (@col) 
	{
		my $results = 0;
		foreach $B (@archiu2) 
		{
			chomp($B);
			my @FILE2 = split(/\t/,$B);
			my $div = $FILE2[96]/$FILE2[94]; # ins
			my $complex = $FILE2[28]/$FILE2[27] if ($FILE2[27] > 0);
			#my $m = $FILE2[3] + $FILE2[46] + $FILE2[69];
			#my $K_wat = $FILE2[95]/$FILE2[93]; 
			##my $watterson = $K_wat/$dom;
			if ($chr eq $FILE2[0] and $FILE2[41] > $recombination1 and $FILE2[41] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $FILE2[100] eq $domain and $complex > $m_complexity1 and $complex <= $m_complexity2)
			{
				my @site_frq = $FILE2[$col_sfs];
				my @site_frq_cell = split(/:/,$site_frq[0]);
				$results += $site_frq_cell[$cell];
			}
		}
		print $results.":";
	}
	print "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0";
}
print "\n";	

exit