#!/usr/bin/perl -w
#  Author: Paul B. Horton, originally adapted from PSORT code by Kenta Nakai
#  Organization: Computational Biology Research Center, AIST, Japan
#  Creation Date: : 2006.4.28
#  Last Modification Date: $Date: 2006/07/11 09:26:19 $
#  Copyright: Paul B. Horton
#  Licensing: Anyone may distribute this under the GNU public license. For other
#             license information contact: horton-p AT aist.go.jp
#
#  Description: Map amino acid residues to integers.
#
#  Purpose: To be used with revised psort code.
#
package PsortAminoAcidIndexer;

use vars qw($VERSION);

$VERSION = "1.00";


my %aminoAcidOneLetterCodeToIndex = (
     'A'=> 0, 'R'=> 1, 'N'=> 2, 'D'=> 3, 'C'=> 4,
     'Q'=> 5, 'E'=> 6, 'G'=> 7, 'H'=> 8, 'I'=> 9,
     'L'=>10, 'K'=>11, 'M'=>12, 'F'=>13, 'P'=>14,
     'S'=>15, 'T'=>16, 'W'=>17, 'Y'=>18, 'V'=>19,
     'B'=>20, 'Z'=>21, 'X'=>22
);


sub toIndex{
    my $aminoAcidOneLetterCode = shift;
    $aminoAcidOneLetterCode || die "Error: toIndex method expects one argument but got none\n";
    if(  defined( $aminoAcidOneLetterCodeToIndex{ $aminoAcidOneLetterCode } )  ){
	$aminoAcidOneLetterCodeToIndex{ $aminoAcidOneLetterCode };
    }else{
	die "Error: unrecognized one letter amino acid code: \"$aminoAcidOneLetterCode\"\n";
    }
}

sub size{
    scalar keys %aminoAcidOneLetterCodeToIndex;
}


1;
    



