#!/usr/bin/perl

#screen_list.pl
#J. Chapman

$n_args = @ARGV;
if (($n_args != 2) && ($n_args != 3)) {
    print "Filter a fasta file for/against a list of entry names.\n";
    print "Usage: ./screen_list.pl <list> <fasta file> <<keep?>>\n";
    print "If a third argument is given the list entries will be
kept,\n";
    print "otherwise they will be rejected.\n";
    exit;
}

my $invert_sense = 0;
if ($n_args == 3) {$invert_sense = 1;}

open(F,$ARGV[0]) || die "Couldn't open file $ARGV[0]\n";
my %bad_ones = ();
while (my $i = <F>) {
    chomp $i;
    my ($name) = $i =~ /(\S+)/;
    $bad_ones{$name} = 1;
}
close F;

open(F,$ARGV[1]) || die "Couldn't open file $ARGV[1]\n";
my $good_one = 1;
while (my $i = <F>) {
    if ($i =~ />/) {
        my ($id) = $i =~ />(\S+)/;
        if (exists($bad_ones{$id})) {$good_one = 0;}
        else {$good_one = 1;}
    }
    if (($invert_sense == 0) &&  ($good_one == 1)) {print $i;}
    elsif (($invert_sense == 1) &&  ($good_one == 0)) {print $i;}
}
close F;
