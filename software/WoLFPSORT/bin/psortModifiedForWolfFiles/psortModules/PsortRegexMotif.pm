#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright: Paul B. Horton 2006, All rights reserved.
#  Creation Date: 2006.8.26
#  Last Modification Date: $Date: 2006/08/28 05:50:05 $
#
#  Description: Motif represented as a regular expression
#
#  Purpose: Originally created which cleaning the PSORT II code.
#
package PsortRegexMotif;
$VERSION = "1.00";
use PsortMotifMatch;
use strict;

use constant DESCRIPTION => 0;
use constant PATTERN => 1;

use overload
    q{""} => 
    sub{
	my $d = $_[0]->description();
	my $p = $_[0]->pattern();
	"description: $d\npattern: $p\n" };

sub new{
    my ($class, $description, $pattern) = @_;
    return bless [$description, $pattern], ref($class) || $class;
}

sub description{
    return $_[0]->[DESCRIPTION];
}

sub pattern{
    return $_[0]->[PATTERN];
}

# non-overlapping matches.
sub matches{
    my ( $self, $seq ) = @_;

    my $pat = $self->pattern();
    my @matchString = $seq =~ /$pat/g;
    my @retVal = ();
    my $searchStartPos = 0;
    for my $matchString (@matchString){
	my $pos = index( $seq, $matchString, $searchStartPos );
	$searchStartPos = $pos + 1;
	push @retVal, PsortMotifMatch->new( $self->description(),
					    $pos, 
					    length($matchString) );
    }
    return @retVal;
}


# non-overlapping matches.
sub numMatches{
    my ( $self, $seq ) = @_;

    my $pattern = $self->pattern();
    my @matches = $seq =~ /$pattern/g;
    return scalar @matches;
}

1;
