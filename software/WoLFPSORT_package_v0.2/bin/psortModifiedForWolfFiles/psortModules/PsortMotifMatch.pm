#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Copyright: Paul B. Horton 2006, All rights reserved.
#  Creation Date: 2006.8.26
#  Last Modification Date: $Date: 2006/08/28 05:46:38 $
#
#  Description: Motif match to the current sequence.
#
#  Indexing: indices start from 0. 
#            intervals are half open as in the C++ STL.
#            i.e. the first base of the motif match is $seq[begIdx()]
#            and the last base is $seq[endIdx()-1]
#
#  Purpose: Originally created which cleaning the PSORT II code.
#
use strict;
package PsortMotifMatch;

use vars qw($VERSION);
$VERSION = "1.00";

use constant DESCRIPTION => 0;
use constant OFFSET => 1;
use constant LENGTH => 2;

sub new{
    my ($class, $description, $begIdx, $endIdx ) = @_;
    my $data = [$description, $begIdx, $endIdx];
    return bless $data, ref($class) || $class;
}

sub description{
    return $_[0]->[DESCRIPTION];
}

sub offset{
    return $_[0]->[OFFSET];
}

sub length{
    return $_[0]->[LENGTH];
}
