#!/usr/bin/perl -w
#
#  Author: "Paul B. Horton"
#  Organization: Computational Biology Research Center, AIST, Japan
#  Creation Date: 2004.12.9
#  Last Modified: $Date: 2006/08/31 08:25:13 $
#  Copyright: Paul Horton
#  License: See license information for WoLF PSORT package
#  
#  Description: Read sequence input from standard in. Check for problems and clean up
#               if necessary or output error warning.
#
#  Purpose: Originally developed for use with WoLF PSORT.
#
#  Input: multifasta format sequences stream, possibly with errors.
#
#  Output: A multifasta format file containing the sequences which
#          could be parsed (possibly upcased or with some other format changes),
#          and also warning or error lines. Warning and error lines start
#          with "#"
#
use strict;
package checkFastaInput;

sub setId();
sub printSeq( $ );
sub checkAndPrintRecord();

# global variables.
my $minSeqLen = 30;
my $maxSeqLen = 10000;
my $seq = "";
my $idLine = "";
my $id = "";
my $recordHeadPattern = '^\s*>'; # allows white space before '>'
my $usage = "$0 [--html]";

my $endl = "\n";  # end of line character to use for diagnostic messages.
my $arg = shift;

if( $arg ){
    if( $arg eq "--html" ){
	$endl = "<BR>\n";
    }else{
	print "Error, unparsed arg: $arg\n\n$usage\n";
	exit -1;
    }
}
    

{ # Main.
    my $line = <STDIN>;
    $line || die "Error: empty input$endl";
    $line =~ tr/\015//d; # get rid of nasty ^M's
    while( $line =~ /^\#/ ){
	$line = <STDIN>;  # skip comment lines.
	$line =~ tr/\015//d; # get rid of nasty ^M's
    }
    # if raw sequence given, add id line. No way to handle multiple raw sequences.
    if( $line =~ /$recordHeadPattern/o ){
	$idLine = canonicalizeRecordHeadLine( $line );
    }else{
	if( $seq =~ />/ ){
	    print "#Error: sequence contains a '>' character, which indicates bad input$endl";
	    print "#Error: seq: \"$seq\"$endl";
	    exit( -1 ); # might as well exit here because the input is so weird.
	}
	$idLine = "> queryProtein";
	$line =~ tr/a-z/A-Z/;
	$line =~ tr/A-Z//cd;
	$seq .= $line;
    }

    setId(); # get first id.

    # loop through entries.
    while( <STDIN> ){
	$line = $_;
	$line =~ tr/\015//d; # get rid of nasty ^M's
	next  if( $line =~ /^\#/ );  # skip comment line.
	if( $line =~ /$recordHeadPattern/o ){
	    checkAndPrintRecord();
	    $seq = "";
	    $idLine = canonicalizeRecordHeadLine( $line );
	    setId();
	}else{
	    $line =~ tr/a-z/A-Z/;
	    $line =~ tr/A-Z//cd;
	    $seq .= $line;
	}
    }
    
    # handle last record
    checkAndPrintRecord();
} # end main.

exit 1; 


{   my $unlabeledQueryCount = 0;
    sub setId(){
	if( $idLine =~ /$recordHeadPattern\s*(\S+)\s*/ ){
	    $id = $1;
	}elsif( $idLine =~ /$recordHeadPattern\s*$/ ){  #  Line starts with '>' but has only white space after that.
	    ++$unlabeledQueryCount;    #  in this case give ids named: queryProtein1, queryProtein2, ...
	    $id = "queryProtein$unlabeledQueryCount";
	    $idLine = ">$id";
	}else{
	    print "#Error: Reaching this point implies a logical error in this program: idLine was:", $idLine, "\"$endl";
	    exit -1;
	}
    }
}


sub canonicalizeRecordHeadLine{
    my $idLine = shift;
    chomp $idLine;
    $idLine =~ s/$recordHeadPattern/>/og;
    return $idLine;
}
    

sub checkAndPrintRecord(){
    if( length($seq) < $minSeqLen ){
	print "#Warning: Sequence $id skipped because length ", length($seq), " less than $minSeqLen$endl";
	return;
    }
    if( length($seq) > $maxSeqLen ){
	print "#Warning: Sequence $id skipped because length ", length($seq), " greater than $maxSeqLen$endl";
	return;
    }
    if( $seq =~ /O/ ){
	print "#Warning: Sequence $id skipped because it contains a non-standard amino acid 'O'$endl";
	print "#seq is: $seq$endl";
	return;
    }
    ( $seq =~ /^M/ ) || print "#Warning: sequence $id does not start with M$endl";
    print "$idLine\n";
    printSeq( $seq );
}

sub printSeq( $ ){
    my $maxSeqLineLen = 60;
    for( my $start = 0; $start < length( $seq ); $start += $maxSeqLineLen ){
	my $curLineLen = length( $seq ) - $start;
	$curLineLen = $maxSeqLineLen  if( $curLineLen > $maxSeqLineLen );
	print substr( $seq, $start, $curLineLen ), "\n";
    }
}
