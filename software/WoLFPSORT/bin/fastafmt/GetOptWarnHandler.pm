#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, Japan
#  Creation Date: : 2006.6.9
#  Last Modification Date: $Date: 2006/08/26 10:21:17 $
#  Copyright: (C) Paul B. Horton 2006. All Rights Reserved.
#
#  Description: A few functions for use with the standard module Getopt::Long
#               Using these functions you can tell the user exactly which option
#               failed to parse correctly.
#
#  Usage: In the program doing the ARGV parsing with Getopt::Long
#         insert the following kind of code:
#
#         fastafmt::GetOptWarnHandler->before( @ARGV );
#         my $getOptionsRetval = GetOptions( ...
#         fastafmt::GetOptWarnHandler->after( @ARGV );
#
use strict;

package fastafmt::GetOptWarnHandler;

my @origARGV;

sub before{
    @origARGV = @_;
    $SIG{__WARN__} = \&getOptWarnHandler;
}

sub after{
    $SIG{__WARN__} = 'DEFAULT';
    my @unparsedArg = grep {/^--\S+|^-\S$/} @ARGV;
    if( @unparsedArg ){
	print "Error: could not parse the following options\n";
	for( my $i = 0; $i < @unparsedArg; ++$i ){
	    $i && print ", ";
	    print "$unparsedArg[$i]";
	}
	print "\n\n";
    }
    return scalar @unparsedArg;			     
}

sub getOptWarnHandler{
    while( @origARGV ){
	my $curOrigArg = pop @origARGV;
	if( !@ARGV || ($curOrigArg ne pop @ARGV) ){
	    print "Error in parsing command line. Perhaps problem with \"$curOrigArg\"\n";
	    last;
	}
    }
    ::pod2usage(0);
    exit( -1 );
}


1;
