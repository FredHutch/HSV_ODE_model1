#!/usr/bin/perl
my $low=`grep c_low hsv_sim.in | cut -d' ' -f2 | tail -1`;
my $high= `grep c_high hsv_sim.in | cut -d' ' -f2 | tail -1`;
my $span;
my $init=`grep 'Total Score' $ARGV[0] | cut -d' ' -f 4,8 | sort -n -r | head -1 | cut -d'=' -f 2`;
my $squeeze=1;
if ($high - $init < ($high - $low)/5.0) {
   $squeeze=0;
}
if ($init - $low < ($high - $low)/5.0) {
   $squeeze=0;
}
if ($squeeze == 1) {
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
printf ("c_init %g\nc_low %g\nc_high %g\n",$init,$low,$high);
