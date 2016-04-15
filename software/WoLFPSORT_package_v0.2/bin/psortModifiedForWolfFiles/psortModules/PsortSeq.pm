#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright (C) 2006, Paul B. Horton, All rights reserved.
#  Creation Date: 2006.8.26
#  Last Modification Date: $Date: 2006/08/27 10:59:19 $
#
#  Description: Global Sequence Variable for PSORT.
#               Note: this implementation is not object oriented.
#               You get exactly one sequence and no encapsulation.
#               I may smoke a turd in hell for that...
#
BEGIN{  
    use FindBin;
    push @INC, "$FindBin::Bin";
}
use PsortFastaReader;
use strict;

package PsortSeq;

# declare package variables
our @seq = ();           # the sequence as a list
our $seq = "";           # the sequence as a string
our $id = "";            # the sequence id
our $loc = "";           # the sequence localization class
our $description = "";   # the sequence description
our $fastaHead = "";     # the head of the sequence's fasta input record
our $len = 0;            # the sequence len

# Whether the the head of the fasta input record is expected to
# have a field after the id representing a localization class
# (this specialization of the fasta format has been called psort format)
our $expectLabeledInput = 0;  

# list of non-standard characters orgininally found in the sequence.
our $badAminoAcidCharsFoundList = ();

# fasta reader.
our $fastaReader = PsortFastaReader->new();


sub readOneSequence{
    my $seqInputFile = shift;

    $fastaReader->getRecord( $seqInputFile ) || return 0;
    $seq = $fastaReader->recordBodyAsOneLine();
    $len = length( $seq );
    $id = $fastaReader->recordId();
    $loc = $expectLabeledInput ? $fastaReader->recordField( 1 ) : "unknown";
    $description = $fastaReader->recordDescription();
    $fastaHead = $fastaReader->recordHead();

    my $weirdCount = ($seq =~ tr/X/G/);
    $weirdCount && warn "$id: treating $weirdCount X's as Glycines\n";

    $weirdCount = ($seq =~ tr/B/D/);
    $weirdCount && warn "$id: treating $weirdCount B's as Aspartate\n";

    $weirdCount = ($seq =~ tr/J/C/);
    $weirdCount && warn "$id: treating $weirdCount J's as Cysteine\n";

    $weirdCount = ($seq =~ tr/U/C/);
    $weirdCount && warn "$id: treating $weirdCount U's as Cysteine\n";

    $badAminoAcidCharsFoundList = "";
    unless( $seq =~ /^[ABRNDCQEGHIJLKMFPSTUWXYV]*$/ ){
	print "# setting seqBadAminoAcidChars, seq:|$seq|\n";
	$badAminoAcidCharsFoundList = $seq;
	$badAminoAcidCharsFoundList =~ tr/ABRNDCQEGHIJLKMFPSTUWXYV//d;
    }

    @seq = split(//, $seq);
} # end readOneSequence.


sub prettyPrint {
    my $charsPerBlock = 10;
    my $numFullBlocks = int( length($seq) / $charsPerBlock );
    if ($numFullBlocks > 0) {
	foreach my $i (1..$numFullBlocks) {
	    print substr( $seq, ($i-1) * $charsPerBlock, $charsPerBlock );
	    print( ( $i % 5 ) ? " " : "\n" );
	}
    }
    print substr( $seq, $numFullBlocks * $charsPerBlock ), "\n";
}


1;
