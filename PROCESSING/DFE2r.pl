#!/usr/bin/perl

use Getopt::Long;
use Cwd;
use warnings;
GetOptions
(
	'file=s'   => \$file,
	'output=s' => \$output,
);
 my $dir = getcwd;
open(INPUT, "$file") or die "Error 1";

@file = ();
@final = ();

while (<INPUT>) 
{
	if ($_ =~ m/^l/g) 
	{
		push(@file, $_);
	}
	else 
	{
		print $_;
	}
}

for (@file) 
{
   s/(lambda -?\d+.\d+ selected_divergence -0.000000 alpha inf omega_A (-?\d+\.\d+))/-0.000000\tNA\t$2\tNA\tNA\t0.000000/g;
   s/(lambda -?nan selected_divergence -?nan alpha -?nan omega_A -?nan)/NA\tNA\tNA\tNA\tNA\tNA/g;
   s/(lambda -?\d+.\d+ selected_divergence )//g;
   s/( alpha | omega_A )/\t/g;
}

@a = ( ['2L','2R','3L','3R','X'], 
       ['low','medium','high'],   
       ['active','both','inactive'],  
       ['low','medium','high','sup'], 
       ['low','high']  
);

$header = "Ka\talpha\tomega_A\tKa_pos\tKa_neg\tKs\n";


print "Feature 1: ";
$feature1 = <>;
chomp($feature1);
#print "Feature 2: ";
#$feature2 = <>;
#chomp($feature2);

$n=1;
$s=1;

if ($feature2 eq "domain") {
	$n=2;
}
if ($feature1 eq "domain")
{
	$s=2;
}
if ($feature1 eq "order")
{
	$s=3;
}
if ($feature2 eq "order")
{
	$n=3;
}
if ($feature1 eq "erep")
{
	$s=4;
}
if($feature2 eq "erep")
{
	$n=4;
}
if ($feature1 eq "order2")
{
	$s=4;
}
if($feature2 eq "order2")
{
	$n=4;
}

foreach $i (@{$a[0]}) {
    foreach $c (@{$a[1]}) {
        foreach $k (@{$a[1]}) {
        	foreach $j (@{$a[$s]}) {
        		#foreach $m (@{$a[$n]}) {
            	$arr = $i."\t".$c."\t".$k."\t".$j."\t";#.$m."\t";
            	push(@final, $arr);
            	#}
        	}
        }
    }
}

for (@file) 
{
	my @vareach=split("\t",$_);
    if ($_ =~ /^-0.000000/ )
    {
    	push(@comb,$_);
    	break;
    }
    elsif ($_ !~ /^NA/ )
    {
		$Ka_pos = $vareach[0] * $vareach[1];
		$Ka_pos = sprintf "%.10f", $Ka_pos;
		$Ka_neg = $vareach[0] - $Ka_pos;
		$Ka_neg = sprintf "%.10f", $Ka_neg;
		$Ks = $vareach[0]/$vareach[2];
		$Ks = sprintf "%.10f", $Ks;
		chomp($_);
		$comb = $_."\t$Ka_pos\t$Ka_neg\t$Ks\n";
		push(@comb, $comb);
	}
	else 
	{
   	push(@comb,$_);
    }
}
open(OUTPUT,">$output") or die "Error 2";
$i=0;
#print OUTPUT "chr\tmutation\trecombination\t$feature1\t$feature2\t$header";
print OUTPUT "chr\tmutation\trecombination\t$feature1\t$header";
print "chr\tmutation\trecombination\t$feature1\t$feature2\t$header\n";
for (@final) {
	$kk = $final[$i].$comb[$i];
	print OUTPUT $kk;
	print $kk;
	$i++;
}
print OUTPUT "\n";
print "\n\nFile $output saved in $dir.\nFeatures evaluated: chr, mutation, recombination, $feature1, $feature2\n";

exit