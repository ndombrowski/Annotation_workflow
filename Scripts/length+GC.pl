#! /usr/bin/perl

#length+GC.pl
#Gene W. Tyson

### Input files ###

$n_args = @ARGV;

if ($n_args != 1) {
    print "\nThis program takes a fasta file, extracts length\n"; 
    print "and %GC information\n";
    print "Usage: ./length+GC.pl\n\n";
    print "e.g ./";
    exit;
}

open (CONTIGS, $ARGV[0]) || die "Couldn't open $ARGV[1]\n";

$/= ">";

while (my $b = <CONTIGS>) {
    chomp $b;
    next unless $b;
    my ($name, @sequence) = split (/\n/, $b);
    my $seq = join ("", @sequence);
    $length = length($seq);
    push (@names, $name);
    $sequences{$name} = uc $seq;
}
close CONTIGS;

foreach $value (@names) {
    $seq = $sequences{$value};
    @sequence = (split //, $seq);
    my $size = @sequence;
    
    foreach my $base (@sequence) {
	if ($base eq "G") {$g++;}
	if ($base eq "C") {$c++;}
    }
    my $GC = (($g+$c)/$size);
    $GC_conversion = (int($GC*1000))/1000;
    print "$value\t$GC_conversion\t$size\n";
    $g=0;
    $c=0;
    $size=0;
}
exit;
