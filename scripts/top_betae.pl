#!/usr/bin/perl
my $low=`grep betae_low hsv_sim.in | cut -d' ' -f2 | tail -1`;
my $high= `grep betae_high hsv_sim.in | cut -d' ' -f2 | tail -1`;
my $span= ($high - $low)/4.0;
my $init=`grep 'Total Score' $ARGV[0] | cut -d' ' -f 4,14 | sort -n -r | head -1 | cut -d'=' -f 2`;
$squeeze=1;
if ($high - $init < ($high - $low)/5.0) {
   $squeeze=0;
}
if ($init - $low < ($high - $low)/5.0) {
   $squeeze=0;
}
if (defined($squeeze) && $squeeze == 1) {
   $span= ($high - $low)/4.0;
} else {
   $span= ($high - $low)/2.0;
}
$low=$init-$span;
if ($low < 0) {
    $low=0.0;
    $span=2.0*$init;
}
$high=$init+$span;
printf ("betae_init %g\nbetae_low %g\nbetae_high %g\n",$init,$low,$high);
