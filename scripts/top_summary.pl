#!/usr/bin/perl
my $top_file;
if ( $ARGV[0] ne "") {
    $top_file=$ARGV[0];
} else {
    die "Usage: top_summary <filename>";
}

if ( ! -e $ARGV[0]) {
    die "The file $ARGV[0] does not exist";
}

my $top_score=`./top_scores.sh $top_file | cut -f 1 -d' '`;
chomp $top_score;
chomp $top_file;
print "Top score $top_score from $top_file";
print "\n";
$result=`grep -C 66 "Total Score = $top_score" $top_file | head -67`;
print $result;
