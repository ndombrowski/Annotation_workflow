#!/usr/bin/perl -w

# Jennifer Meneghin
# 08/14/2007
# This script removes duplicate records from a "short" format BLAST output file, and keeps only the "best" records  (sorts by smallest e-value and then biggest percent identity)
# Usage: blast_best.pl <input file> <output file>

#-----------------------------------------------------------------------------------------------------------------------------------------
#Deal with passed parameters
#-----------------------------------------------------------------------------------------------------------------------------------------
#If no arguments are passed, show usage message and exit program.
if ($#ARGV == -1) {
    usage("BLAST BEST 1.0 2007");
    exit;
}

#get the names of the input file (first argument passed) and output file (second argument passed)
$in_file = $ARGV[0];
$out_file = $ARGV[1];

#Open the input file for reading, open the output file for writing.
#If either are unsuccessful, print an error message and exit program.
unless ( open(IN, "$in_file") ) {
    usage("Got a bad input file: $in_file");
    exit;
}
unless ( open(OUT, ">$out_file") ) {
    usage("Got a bad output file: $out_file");
    exit;
}

#Everything looks good. Print the parameters we've found.
print "Parameters:\ninput file = $in_file\noutput file = $out_file\n\n";

#-----------------------------------------------------------------------------------------------------------------------------------------
#The main event
#-----------------------------------------------------------------------------------------------------------------------------------------

$counter = 0;
$total_counter = 0;

print "De-duplicating File...\n";

@in = <IN>;

#Do stuff for each line of text in the input file.
foreach $line (@in) {
    #if the line starts with a pound symbol, it is not real data, so skip this line.
    if ( $line =~ /^#/ ) {
	 next;
     }

    #Count the total number of data lines in the file.
    $total_counter++;

    #The chomp commands removes any new line (and carriage return) characters from the end of the line.
    chomp($line);

    #Split up the tab delimited line, naming only the variables we are interested in.
    ($id, $junk,  $junk, $junk, $evalue, $bitscore, $junk, $junk, $junk, $junk, $junk, $junk, $junk) = split(/\t/, $line);

    #check to see if the id label is already in the list of ids (called dedupe)
    #if its not there, add it.
    if ( $dedupe{$id} ) {
	#if it is, look at the old line to see if it is still "better" than the new one.
	( $junk, $junk,  $junk, $junk, $list_evalue, $list_bitscore, $junk, $junk, $junk, $junk, $junk, $junk, $junk) = split(/\t/,$dedupe{$id});

	#if the new evalue is better than the old one, change the value of this id to the new line.
	#otherwise, if the the new evalue is the same, and the percent_identity is better, change the value of this id to the new line.
	#otherwise, don't do anything (keep the old line).
	if ( $evalue < $list_evalue ) {
	    $dedupe{$id} = $line;
	}
	elsif ( $evalue == $list_evalue ) {
	    if ( $bitscore > $list_bitscore ) {
		$dedupe{$id} = $line;
	    }
	}
    }
    else {
	$dedupe{$id} = $line;
	#count the number of non-duplicated lines we have.
	$counter++;
    }
}
print "Total # records = $total_counter\nBest only # records = $counter\n";
print "Writing to output file...\n";

#Print the final "dedupe" list to the new file (adding the new line back on the end).
foreach $id (sort keys %dedupe) {
    print OUT "$dedupe{$id}\n";
}

#Close the files.
close(IN);
close(OUT);
print "Done.\n";

#-----------------------------------------------------------------------------------------------------------------------------------------
#Subroutines
#-----------------------------------------------------------------------------------------------------------------------------------------
sub usage {
    my($message) = @_;
    print "\n$message\n";

    print "\nThis script removes duplicate records from a \"short\" format BLAST output file, and keeps only the \"best\" records.\nIt sorts by smallest e-value and then biggest percent identity.\n";
    print "Usage: blast_best.pl <input file> <output file>\n";
    print "\nJennifer Meneghin\n";
    print "08/15/2007\n";
}

