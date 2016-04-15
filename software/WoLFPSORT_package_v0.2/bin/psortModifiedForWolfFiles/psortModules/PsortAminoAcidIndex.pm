#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Creation Date: : 2006.4.28
#  Last Modification Date: $Date: 2006/07/11 09:14:26 $
#  Copyright: Paul B. Horton
#  Licensing: Anyone may distribute this under the GNU public license. For other
#             license information contact: horton-p AT aist.go.jp
#
#  Description: Represent a amino acid index.
#               Basically a mapping from amino acids to real numbers.
#
#               This class has many methods for calculating the average
#               or maximum value of a fixed length window within a
#               sequence. From a design standpoint it would be cleaner
#               if that functionality were factored out into another class
#               or perhaps a specialization of this class.
#
# Notes: In this class, sequences are represented as a reference to an
#        array which has one residue per element. The variable name
#        used for this is $seqArrayRef;
#
#  Purpose: Created for use with WoLF PSORT.
#
use strict;
package PsortAminoAcidIndex;

use vars qw($VERSION);

$VERSION = "1.00";


########################################################################
#
# Amino Acid Indices Data
#
# B was calculated as an abundance weighted average of aspartic acid (D) and asparagine (N)
# Z was calculated as an abundance weighted average of glutamine (Q)  and glutamic acid (E)
# X was calculated as an abundance weighted average of the 20 standard amino acids.
# U for selenocysteine is not currently included because PSORT treats it as cysteine.
#
# The relative abundances used are given as the first index.
#
# from NCBI Amino Acid Explorer web site 2006/7/11
my %relativeAbundance =
    ( name => "relative abundance",
      L	=> 0.0994, A => 0.0884, G => 0.0703, V => 0.0677, S => 0.0672,
      E	=> 0.0624, I => 0.0595, R => 0.0570, T => 0.0543, D => 0.0539,
      K	=> 0.0527, P => 0.0471, N => 0.0417, F => 0.0400, Q => 0.0382,
      Y => 0.0300, M => 0.0237, H => 0.0220, C => 0.0124, W => 0.0121 );

my %EngelmanGES1986 = 
    ( name => "Engelman GES 1986",
      A =>  -6.7, R =>  51.5, N => 20.1, D =>  38.5, C =>  -8.4, Q => 17.2,
      E =>  34.3, G =>  -4.2, H => 12.6, I => -13.0, L => -11.7, K => 36.8,
      M => -14.2, F => -15.5, P =>  0.8, S =>  -2.5, T =>  -5.0, W => -7.9,
      Y =>   2.9, V => -10.9, B => 30.5, Z =>  27.8, X =>   5.8 );

# was called "hydro" in old PSORT code
my %KYTJ820101 =
    ( name => "Kyte Doolittle",
      id   => "KYTJ820101",
      A =>  1.8, R => -4.5, N => -3.5, D => -3.5, C =>  2.5, Q => -3.5,
      E => -3.5, G => -0.4, H => -3.2, I =>  4.5, L =>  3.8, K => -3.9,
      M =>  1.9, F =>  2.8, P => -1.6, S => -0.8, T => -0.7, W => -0.9,
      Y => -1.3, V =>  4.2, B => -3.5, Z => -3.5, X => -0.19 );

my %KLEP840101 =
    ( name => "Klein net Charge",
      id   => "KLEP840101",
      A =>  0, R =>  1, N => 0,     D => -1,     C => 0, Q => 0,
      E => -1, G =>  0, H => 0,     I =>  0,     L => 0, K =>  1,
      M =>  0, F =>  0, P => 0,     S =>  0,     T => 0, W => 0,
      Y =>  0, V =>  0, B => -0.56, Z =>  -0.62, X => -0.01 );


my %nameToReference = 
    ( engelmanges1986 => \%EngelmanGES1986,
      kytedoolittle   => \%KYTJ820101,
      kleinnetcharge  => \%KLEP840101 );


my %idToReference = 
    ( engelmanges1986 => \%EngelmanGES1986,
      KLEP840101      => \%KLEP840101,
      KYTJ820101      => \%KYTJ820101 );


####################################################################
#
# Object Methods
#
sub new{
    my $indexIdOrName = $_[1];
    if( _idToReference($indexIdOrName) ){
	bless _idToReference( $indexIdOrName ), $_[0];
    }elsif( _nameToReference($indexIdOrName) ){
	bless _nameToReference( $indexIdOrName ), $_[0];
    }else{
	die "No amino acid with id or name: \"$indexIdOrName\" known\n";
    }
}

sub name{  $_[0]{ "name" };  }
sub id  {  $_[0]{ "id" };    }

sub value{
    ( @_ == 2 ) || die "expected 2 arguments\n";
    my( $this, $residue ) = @_;
   
    (length( $residue ) == 1) || die "\"$residue\" not a valid amino acid residue\n";
    exists( $$this{$residue} ) || die "\"$residue\" not a valid amino acid residue\n";
    return $$this{$residue};
}


my @standard20Residues = qw( A C D E F G H I K L M N P Q R S T V W Y );


sub abundanceWeightedAverageValue{
    ( @_ == 1 ) || die "expected one argument";
    my $this = shift;

    my $sum = 0;

    my $abundanceSum = 0;
    for my $res (@standard20Residues){
	exists( $relativeAbundance{$res} ) || die "no abundance data for res: \"$res\"\n";
	$sum += $relativeAbundance{$res} * $this->value($res);
	$abundanceSum += $relativeAbundance{$res};
    }
    return( $sum / $abundanceSum );
}


######################################################################
#
# Sequence or sequence window scanning methods which perhaps belong
# more in a specialization that in a base amino acid index class.
#
######################################################################

# average value for residues from a seq array.
sub seqAverage{
    my $this = shift;
    my $seqArrayRef = $_[0];
    my( $startPos, $endPos ) = _startEndPos( @_ );
    
    my $sum = 0;
    for( my $i = $startPos; $i < $endPos; ++$i ){
	$sum += $this->value( $$seqArrayRef[$i] );
    }
    return $sum;
}


# return the maximum of the average index value of all $windowLength length substrings in
# the range given by $offset and $length
sub maxWindowAverage{
    my $this = shift;
    my $seqArrayRef = shift;
    my $windowLength = shift;
    unshift @_, $seqArrayRef;
    my( $startPos, $endPos ) = _startEndPos( @_ );

    my( $numWindowPositions ) = $endPos - $startPos - $windowLength + 1;
    if( $numWindowPositions < 0 ){
	die "startPos: $startPos, endPos: $endPos, windowLength: $windowLength";
    }
    my $max = $this->seqAverage( $seqArrayRef, $startPos, $windowLength );
    my( $average ) = $max;
    for( my $i = $startPos; $i < $numWindowPositions; ++$i ){
	$average -= $this->value( $seqArrayRef->[$i] );
	$average += $this->value( $seqArrayRef->[$i+$windowLength] );
	$max = $average if( $average > $max );
    }
    return $max;
}



# translate substr or splice style seq, offset, length
# arguments into starting and ending positions.
sub _startEndPos{
    @_ || die "expected arguments";
    my $seqArrayRef = shift;
    (@_ == 0) && push @_, 0;
    (@_ == 1) && push( @_, scalar @{$seqArrayRef} );
    (@_ == 2) || die "expected 1 to 3 arguments";
    my( $startPos, $length ) = @_;
    my $seqLen = scalar @{$seqArrayRef};

    if( $startPos < 0 ){       die "invalid startPos value"; }
    if( $startPos > $seqLen ){ die "invalid startPos value"; }
    my $endPos = $startPos + $length;
    if( $length < 0 ){
	$endPos = $seqLen - $length;
    }
    if( $endPos < $endPos ){   die "invalid endPos value"; }
    if( $endPos > $seqLen ){   die "endPos: $endPos, > seqLen: $seqLen"; }
    return( $startPos, $endPos );
}

    



########################################
#
# Public class methods
#
sub indexNames{
    return( keys %nameToReference );
}


########################################
#
# Private class methods
#
sub _nameToReference{
    my $indexName = shift;
    $indexName =~ tr/ //d;
    $indexName =~ tr [A-Z] [a-z];
    if(  exists( $nameToReference{$indexName} )  ){
	$nameToReference{$indexName}
    }else{
	return 0;
    }
}

sub _idToReference{
    my $id = shift;
    if(  exists( $idToReference{$id} )  ){
	$idToReference{$id};
    }else{
	return 0;
    }
}

sub _printSeqArrayRef{
    my $seqArrayRef = shift;
    my $ostream = shift;
    $ostream || ( $ostream = *STDOUT);
    
    print $ostream join( "", @{$seqArrayRef} );
}

1;
