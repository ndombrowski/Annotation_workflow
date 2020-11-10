use strict;
use warnings;

my $treeFile = pop;
my %taxonomy = map { /(\S+)\s+(.+)/; $1 => $2 } <>;

push @ARGV, $treeFile;

while ( my $line = <> ) {
    $line =~ s/\b$_\b/$taxonomy{$_}/g for keys %taxonomy;
    print $line;
}
