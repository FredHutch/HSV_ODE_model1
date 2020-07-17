#!/usr/bin/perl
my $line;
my $runnum=0;
my $crit1bin=0;
my $crit1bins=8;
my $crit11bin=0;
my $crit11bins=40;

my $outfile="all.csv";

my %cmax;
my %gamma;
my %betaun;
my %m;
my %ic50;
my %crit11;
my %crit1;
my %episodes;
my %transmissions;
my %acts;
my %positive;
my %over1mm;
my %over2mm;

open(OUTF,">$outfile") || die " Could not open output file $outfile!\n";

if (!defined($ARGV[0])) {
    die "This script requires atleast 1 input file\n";
}

$transmissions[0]=0;
$acts[0]=0;
for (my $i=0; $i < $#ARGV+1; $i++) {

    open(INF,"<$ARGV[$i]") || die " Could not open input file $ARGV[$i]!\n";

    while( $line=<INF>) {
	chomp($line);              # remove the newline from $line.

	my @pieces = split(/ /,$line);

	# Total Score = 0.224477 beta=5.400000e-8 vburstrate=82.000000 p=102991.171934 c=8.800000 theta=2.840000 delta=0.001470 r=42.000000 inf=194000.000000 rinf=295000.000000 beta_e=2.150000e-11 rho=0.690000 eclipse=0.960000 beta_un=0.879000e-11 absorb=0.125000 gamma=5.540000 Cmax=1.400000 IC50=0.210000 m=2.000000 time=106.044000
	if ($pieces[0] =~ /Total/ && $pieces[1] =~ /Score/ ) {
		my @subpieces = split(/=/,$pieces[19]); # cmax
		$cmax[$runnum]=$subpieces[1];

		@subpieces = split(/=/,$pieces[18]); # gamma
		$gamma[$runnum]=$subpieces[1];

		@subpieces = split(/=/,$pieces[20]); # ic50
		$ic50[$runnum]=$subpieces[1];

		@subpieces = split(/=/,$pieces[21]); # m
		$m[$runnum]=$subpieces[1];

		@subpieces = split(/=/,$pieces[16]); # betaun
		$betaun[$runnum]=$subpieces[1];

	} elsif ($pieces[0] eq "criteria" && $pieces[1] =~ /11/) {
#	criteria 11[1]: 4.000000-4.100000, 11719 swabs = 0.321069%
#
		$crit11[$runnum][$crit11bin]=$pieces[6];
		$crit11bin++;

	} elsif ($pieces[0] eq "criteria" && $pieces[1] =~ /1\[/) {
#	criteria 1[2]: 2.884615 vs. mean of 3.140000 - err = 0.255385 * 0.048921
#
		$crit1[$runnum][$crit1bin]=$pieces[2];
		$crit1bin++;

	} elsif ($pieces[0] eq "Percent" && $pieces[5] =~ /1mm/) {
#	Percent Time in episodes with >1mm plaques= 11.783753
		$over1mm[$runnum]=$pieces[7];

	} elsif ($pieces[0] eq "Percent" && $pieces[5] =~ /2mm/) {
#	Percent Time in episodes with >2mm plaques= 2.598658
		$over2mm[$runnum]=$pieces[7];
	} elsif ($pieces[0] eq "No" && $pieces[1] eq "transmissions") {
		$acts[$runnum]=500;
	} elsif ($pieces[0] eq "Transmission" && $pieces[1] eq "after") {
		$transmissions[$runnum]=1;
		$acts[$runnum]=$pieces[2];
	} elsif ($pieces[1] eq "episodes" && $pieces[5] =~ /patient/) {
#3 episodes (tot = 3)for patient 1 (cont epis=2)
		$episodes[$runnum]=$pieces[0];
	} elsif ($pieces[0] eq "Avg" && $pieces[1] =~ /percentage/) {
#Avg percentage of swabs that were positive= 3.773585
		$positive[$runnum]=$pieces[7];

		$crit1bin=0;
		$crit11bin=0;
		$runnum++;
		$transmissions[$runnum]=0;
		$acts[$runnum]=0;
	}

    }
}
close (INF);

printf OUTF ("run,betaun,gamma,ic50,m,cmax,episodes,positive,over 1mm,over 2mm,acts,transmission");
for (my $i=0; $i < $crit1bins; $i++) {
    printf OUTF (",1[%d]",$i+2);
}
printf OUTF ("\n");
my %crit1tot;
my $epitot=0;
my $transtot=0;
my $acttot=0;
my $postot=0;

for (my $j=0; $j < $crit1bins; $j++) {
    $crit1tot[$j] = 0;
}
for (my $i=0; $i < $runnum; $i++) {
    $epitot+=$episodes[$i];
    $transtot+=$transmissions[$i];
    $acttot+=$acts[$i];
    $postot+=$positive[$i];
    printf OUTF ("%d,", $i+1);
    printf OUTF ("%g,", $betaun[$i]);
    printf OUTF ("%g,", $gamma[$i]);
    printf OUTF ("%g,", $ic50[$i]);
    printf OUTF ("%g,", $m[$i]);
    printf OUTF ("%g,", $cmax[$i]);
    printf OUTF ("%g,", $episodes[$i]);
    printf OUTF ("%g,", $positive[$i]);
    printf OUTF ("%g,", $over1mm[$i]);
    printf OUTF ("%g,", $over2mm[$i]);
    printf OUTF ("%d,", $acts[$i]);
    printf OUTF ("%d", $transmissions[$i]);
    for (my $j=0; $j < $crit1bins; $j++) {
	printf OUTF (",%g",$crit1[$i][$j]);
	$crit1tot[$j] += $crit1[$i][$j];
    }
    printf OUTF ("\n");
}
printf OUTF ("\n");
printf OUTF ("Summary\n");
printf OUTF ("%d runs\n",$runnum);
printf OUTF ("%d acts\n",$acttot);
printf OUTF ("%d episodes\n",$epitot);
printf OUTF ("%d transmissions\n",$transtot);
if ($transtot > 0) {
    printf OUTF ("%g acts per transmission\n",$acttot/$transtot);
} else {
    printf OUTF ("No transmissions in %d acts!\n",$acttot);
}
printf OUTF ("%g percent pos swabs\n",$postot/$runnum);
for (my $i=0; $i < $crit1bins; $i++) {
    printf OUTF ("criteria 1[%d] = %g\n",$i+2,$crit1tot[$i]/$runnum);
}
printf OUTF ("\n");
close(OUTF);
