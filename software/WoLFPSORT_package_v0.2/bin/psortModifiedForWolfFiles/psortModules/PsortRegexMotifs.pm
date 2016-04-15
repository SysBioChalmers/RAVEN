#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright: Paul B. Horton 2006, All rights reserved.
#  Creation Date: 2006.8.26
#  Last Modification Date: $Date: 2006/08/28 05:50:44 $
#
#  Description: Collection of Motifs represented as regular expressions
#
#  Purpose: Originally created which cleaning the PSORT II code.
#
package PsortRegexMotifs;
$VERSION = "1.00";
use PsortRegexMotif;
use strict;

use constant NAME => 0;
use constant MOTIFSREF => 1;


sub new{
    my( $class, $istream ) = @_;
    my @data = readMotifsFromTextStream( $istream );
    bless \@data, ref($class) || $class;
}


# accessors
sub name{
    return $_[0]->[NAME];
}

sub motifsRef{
    return $_[0]->[MOTIFSREF];
}


# class method
sub readMotifsFromTextStream{
    my $istream = shift;

    my @motif = ();
    my $description;
    my $pattern;
    my $numLinesRead = 0;
    my $curLine;
    my $name = "";
    while( <$istream> ){
	chomp;
	$curLine = $_;
	next if( $curLine =~ /^\#/ );
	next unless( $curLine =~ /\S+/ );
	
	# if first line, set name
	unless( $name ){
	    $name = $_;
	    next;
	}

	# else take turns reading description and pattern lines.
	if( $numLinesRead == 0 ){
	    $description = $_;
	}else{
	    $pattern = $_;
	    push @motif, PsortRegexMotif->new( $description, $pattern );
	}
	$numLinesRead = ( ++$numLinesRead % 2 );
    }
    if( $numLinesRead ){
	die "could not parse line: $curLine";
    }
    return( $name, \@motif );
}

# returns number of matches reported.
sub reportMatches{
    my ($self, $seq, $endl ) = @_;
    $endl = "\n" unless $endl;

    my $numMatches = 0;
    my $motsRef = $self->motifsRef();
    my $name = $self->name();
    print "checking ", scalar @$motsRef, " $name: ";
    for my $motif (@$motsRef){
	my @matches = $motif->matches($seq);
	@matches || next;
	$numMatches += scalar @matches;

	my $motifDescription = $motif->description();
	print "$endl      $motifDescription:  *** found ***$endl";
	for my $match (@matches){
	    my $matchString = substr( $seq, $match->offset(), $match->length() );
	    printf "         %s at %d$endl", $matchString, $match->offset() + 1; # output counts from 1.
	}
    }
    $numMatches || print " none$endl";
    print "$endl";
    return $numMatches;
}

sub numMatches{
    my ($self, $seq) = @_;

    my $numMatches = 0;

    my $motsRef = $self->motifsRef();
    for my $motif (@$motsRef){
	$numMatches += $motif->numMatches($seq);
    }
    return $numMatches;
}


sub dump{
    my $self = shift;

    print "collection name: ", $self->name(), "\n";
    my $motsRef = $self->motifsRef();
    for my $motif (@$motsRef){
	print "$motif\n";
    }
}

1;
