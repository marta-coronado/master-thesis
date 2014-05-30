#!/usr/bin/perl

use Getopt::Long;

GetOptions(
	'exons=s' => \$file1,
	'shortintrons=s' => \$file2,
	'filter=s' => \$file3,
);

open(INPUT1, "$file1") or die "Error 1";
@archiu1 = <INPUT1>;	
shift(@archiu1);

open(INPUT2, "$file2") or die "Error 2";
@archiu2 = <INPUT2>;
shift(@archiu2);

# open(INPUT3, "$file3") or die "Error 3";
# @archiu3 = <INPUT3>;

# $file3 =~ s/(:)/ /g;
# # print $file3."\n";
# $filters = $file3;

for ($cnt = 63; $cnt >= 0; $cnt--) # to sum the SFS
{
	unshift (@col,$cnt);
	@col = sort {$a <=> $b} @col;
}
	
my @sum0 = ("6","7","8","12"); #for divergence.txt 0fold 				   mdmel_0f	 mdyak_0f	 seg_0f		 div_0f	
my @sum2 = ("83","84","85","89"); #for divergence.txt 2fold                mdmel_2f    mdyak_2f    seg_2f      div_2f   
my @sum4 = ("55","56","57","61"); #for divergence.txt 4fold                mdmel_4f    mdyak_4f    seg_4f      div_4f   
my @sumi = ("111","112","113","114"); #for divergence.txt short introns    mdmel_ins   mdyak_ins   seg_ins     div_ins

my @sfs0 = ("11"); #for sfs.txt 0fold          sfs_0f
my @sfs2 = ("88"); #for sfs.txt 2fold          sfs_2f
my @sfs4 = ("60");  #for sfs.txt 4fold           sfs_4f
my @sfsi = ("115");#for sfs.txt short introns sfs_ins


### FILTERS ###

chomp $file3;
@FILE1 = split(/:/,$file3);

$chr = $FILE1[0];            # chr
$mutation1 = $FILE1[1];      # mutation
$mutation2 = $FILE1[2];      # mutation
$recombination1 = $FILE1[3]; # recombination
$recombination2 = $FILE1[4]; # recombination 

$breadth1 = $FILE1[5];  
$breadth2 = $FILE1[6];

$distance1 = $FILE1[7];  
$distance2 = $FILE1[8];

#$m1 = $FILE1[7];  #size1
#$m2 = $FILE1[8]; #size2

#$expression1 = $FILE1[5];  
#$expression2 = $FILE1[6];

#$breadth1 = $FILE1[5];  
#$breadth2 = $FILE1[6];

#$m1 = $FILE1[5];  #size1
#$m2 = $FILE1[6]; #size2

#$domain = $FILE1[7];   #domain


# $distance1 = $FILE1[7];  
# $distance2 = $FILE1[8];  

# $order1 = $FILE1[5];    # order
# $order2 = $FILE1[6];    # order

# $Er1 = $FILE1[5];     # inclusion levels
# $Er2 = $FILE1[6];     # inclusion levels 

# $mRNA1 = $FILE1[5];   # transcripts 
# $mRNA2 = $FILE1[6];   # transcripts 

# $domain = $FILE1[5];   #domain

# $num1 = $FILE1[7];      #num exons
# $num2 = $FILE1[8];    # num exons

=xx

my $m = $FILE2[6] + $FILE2[83] + $FILE2[55]; #ojo..antes era dif!
and $m > $m1 and $m <= $m2 

my $Exon_Rep = $FILE2[5]/$FILE2[118] if ($FILE2[118] > 0);
and $Exon_Rep > $Er1 and $Exon_Rep <= $Er2 

and $FILE2[5] > $mRNA1 and $FILE2[5] <= $mRNA2

and $FILE2[46] > $order1 and $FILE2[46] <= $order2 

and $FILE2[118] eq $domain 

and $FILE2[34] > $num1 and $FILE2[34] <= $num2

and $FILE2[36] > $distance1 and $FILE2[36] <= $distance2

and $FILE2[37] > $expression1 and $FILE2[37] <= $expression2 

and $FILE2[38] > $breadth1 and $FILE2[38] <= $breadth2 

=cut

### 0-FOLD ###

foreach $col_sum (@sum0) 
{
	my $results = 0;
	foreach $B (@archiu1) 
	{
		chomp($B);
		my @FILE2 = split(/\t/,$B);
		#my $m = $FILE2[6] + $FILE2[83] + $FILE2[55];
		my $div = $FILE2[114] / $FILE2[112]; # ins
		my $div0 = $FILE2[12]; # div en 0 fold 
		#my $Exon_Rep = ($FILE2[5] / $FILE2[35]) if ($FILE2[35] > 0);                                                                                                            																									                                                                                                            
		$results += $FILE2[$col_sum] if ($FILE2[0] eq $chr and $FILE2[49] > $recombination1 and $FILE2[49] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $div0 > 0 and $FILE2[38] > $breadth1 and $FILE2[38] <= $breadth2 and $FILE2[36] > $distance1 and $FILE2[36] <= $distance2); 
	}
	print $results." ";
}

foreach $col_sfs (@sfs0) 
{
	my @site_frq = "";
	my @site_frq_cell = "";
	foreach $cell (@col) 
	{
		my $results = 0;
		foreach $B (@archiu1) 
		{
			chomp($B);
			my @FILE2 = split(/\t/,$B);
			#my $m = $FILE2[6] + $FILE2[83] + $FILE2[55];
			my $div = $FILE2[114] / $FILE2[112]; # ins
			my $div0 = $FILE2[12]; # div en 0 fold 
			#my $Exon_Rep = $FILE2[5]/$FILE2[35] if ($FILE2[35] > 0);
			if ($FILE2[0] eq $chr and $FILE2[49] > $recombination1 and $FILE2[49] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $div <= $mutation2 and $div0 > 0 and $FILE2[38] > $breadth1 and $FILE2[38] <= $breadth2 and $FILE2[36] > $distance1 and $FILE2[36] <= $distance2) 
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

### SHORT INTRONS ###

foreach $col_sum (@sumi) 
{
	my $results = 0;
	foreach $B (@archiu2) 
	{
		chomp($B);
		my @FILE2 = split(/\t/,$B);
		#my $m = $FILE2[6] + $FILE2[83] + $FILE2[55];
		my $div = $FILE2[114] / $FILE2[112]; # ins	
		my $div0 = $FILE2[12]; # div en 0 fold 
		#my $Exon_Rep = $FILE2[5]/$FILE2[35] if ($FILE2[35] > 0);	
		$results += $FILE2[$col_sum] if ($FILE2[0] eq $chr and $FILE2[49] > $recombination1 and $FILE2[49] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $div0 > 0 and $FILE2[38] > $breadth1 and $FILE2[38] <= $breadth2 and $FILE2[36] > $distance1 and $FILE2[36] <= $distance2); 
	}
	print $results." ";
}

foreach $col_sfs (@sfsi) 
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
			#my $m = $FILE2[6] + $FILE2[83] + $FILE2[55];
			my $div = $FILE2[114] / $FILE2[112]; # ins
			my $div0 = $FILE2[12]; # div en 0 fold 
			#my $Exon_Rep = $FILE2[5]/$FILE2[35] if ($FILE2[35] > 0);
			if ($FILE2[0] eq $chr and $FILE2[49] > $recombination1 and $FILE2[49] <= $recombination2 and $div > $mutation1 and $div <= $mutation2 and $div <= $mutation2 and $div0 > 0 and $FILE2[38] > $breadth1 and $FILE2[38] <= $breadth2 and $FILE2[36] > $distance1 and $FILE2[36] <= $distance2) 
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