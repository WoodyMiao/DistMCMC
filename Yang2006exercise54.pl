#!/usr/bin/perl
use strict;
use warnings;

my $n = 948;
my $x = 90;
my $w = 0.1;
my $u = 0.02; # the mu parameter of the exponential prior
my @range = (0, 1); # range of distance
my $d = 0.2; # initial state of distance
#$d = ($range[1]-$range[0])*(rand) + $range[0]; # random initial state
my $noItera = 100000; # number of iterations
my $samFreq = 10; # sampling frequency

my $maxlikeD = -0.75*log(1-4/3*$x/$n);
printf "\nThe p-distance is  %.4f.\n", $x/$n;
printf "The maximum likelihood extimation of the JC69 distance is %.4f.\n", $maxlikeD;
print "\nMCMC is running...\n\n";

open O, ">", "Yang2006exercise54.log";
print O "Sample\tposterior\tlikelihood\tprior\tdistance\n";

my $accept; # number of accepted states
foreach (1 .. $noItera) {
	my $d_old = $d; # old state
	my $d_new = &NewSta($d_old, $w, @range); # new state
	my $likeOld = ($x)*log(0.75-0.75*exp(-4/3*$d_old)) + ($n-$x)*log(0.25+0.75*exp(-4/3*$d_old)); # log likelihood of old state
	my $likeNew = ($x)*log(0.75-0.75*exp(-4/3*$d_new)) + ($n-$x)*log(0.25+0.75*exp(-4/3*$d_new)); # log likelihood of new state
	my $priorOld = log(1/$u) - $d_old/$u; # log prior of old state
	my $priorNew = log(1/$u) - $d_new/$u; # log prior of old state
	my $postOld = $priorOld + $likeOld;
	my $postNew = $priorNew + $likeNew;
	my $accRat; # Acceptance ratio
	if ($postNew > $postOld) {
		$accRat = 1;
	} else {
		$accRat = exp($postNew-$postOld);
	}
	my $r2 = rand;
	if ($r2 < $accRat) {
		$d = $d_new;
		++$accept;
	}
	my $p = $_/$samFreq;
	if ($p == int($p)) {
		if ($d == $d_new) {
			print O "$_\t$postNew\t$likeNew\t$priorNew\t$d\n";
		} else {
			print O "$_\t$postOld\t$likeOld\t$priorOld\t$d\n";
		}
	}
}
close O;

my $accPro = $accept/$noItera; # acceptance proportion 
printf "MCMC done! The acceptance proportion is %.2f%%\n\n", $accPro*100;


sub NewSta {
	my ($OldStat, $SlidWin, $LoBound, $UpBound) = @_;
	my $r1 = rand;
	my $TemStat = $OldStat + ($r1-0.5)*$SlidWin; # Temp state
	my $NewStat; # New state
	if ($TemStat < $LoBound) {
		$NewStat = 2*$LoBound - $TemStat;
	} elsif ($TemStat > $UpBound) {
		$NewStat = 2*$UpBound - $TemStat;
	} else {
		$NewStat = $TemStat;
	}
	return $NewStat;
}
