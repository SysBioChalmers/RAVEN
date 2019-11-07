#
# Date of last CVS update: $Date: 2006/08/27 11:01:28 $
#
############################################################
##  PSORT II:  subroutine package 1 (970725)
##                      copyright: Kenta Nakai 1996, 1997
############################################################

## Sequence Numbering convention:
##   argument: 1, 2, 3,...   internal: 0, 1, 2,...

# This version is for PSORT+WoLF, slightly modified from original May 2004.
# The variable $endl was added by Paul Horton to allow the end of line string
# to be changed (e.g. to "<br>\n").

# 2006/07/12 "our" declaration added for global variables.

# 2006/01/06 Calculation of mtop recoded by Paul Horton in response to bug
#            found by Mr. Katoh. Old variable names were also replaced with
#            more descriptive ones (but just in that section of code).

# 2006/05/02 Amino acid index code was moved into a separate package, also many
#            my declarations were added to allow "use strict"
BEGIN{  
    use FindBin;
    push @INC, "$FindBin::Bin/psortModules";
}
use PsortAminoAcidIndex;
our ( $OPT_V, $endl );
use strict;
package sub1;

## definitions ###################################################

my $kyteDoolittle = PsortAminoAcidIndex->new( "Kyte Doolittle" );
my $ges = PsortAminoAcidIndex->new( "Engelman GES 1986" );

my $MaxPositionOf_CR_End = 10;
my %tms; # not sure but I think mtype_assign may rely on alom2 setting this variable before being called?? PH. 2006/5/1.
my $mostn; # mytpe_assign seems to rely on alom2 setting this variable before being called. PH. 2006/5/1.
my $LSEG = 17; # used to be assigned in alom2 but used elsewhere as well.

## new McGeoch's method ###################################################
sub psg {
#    my $CRMAX = 11;  ## maximum position of CR_end (external numbering )
    my($SEGLEN) = 8;
    my($PSGTHR) = 4.4;

    my($cr_leng, $poschg, $negchg, $ur_leng) = &regions;

    my $ur_leng2 = ($cr_leng+$ur_leng > 30) ? 30 - $cr_leng : $ur_leng;
    return -10000 if $cr_leng+$SEGLEN > $PsortSeq::len;

    my($ur_peak) = &hydseg2($SEGLEN, $cr_leng+1, $ur_leng2);

    if($OPT_V) {
	print "PSG:  a new signal peptide prediction method$endl";
    }
## ad hoc rules ##
    if($ur_leng >= 60 || ($poschg-$negchg) < 0) {
	$ur_peak = 0;
	print ">>> too long uncharged segment: $ur_leng$endl"
	    if($ur_leng >= 60 && $OPT_V);
	printf ">>> negative net charge in the N-region: %3.0f$endl", 
	$poschg-$negchg if($poschg-$negchg < 0 && $OPT_V);
    }
    if($OPT_V) {
	printf("      N-region:  length %d;  pos.chg %d;  neg.chg %d$endl",
	       $cr_leng, $poschg, $negchg);
	printf("      H-region:  length %d;  peak value %6.2f$endl",
	       $ur_leng, $ur_peak);
	printf("      PSG score: %6.2f$endl$endl", $ur_peak - $PSGTHR);
    }
	
    $ur_peak - $PSGTHR;
}

# specifically for gram-positive bacteria
sub psggp {
#	my $CRMAX = 11;  ## maximum position of CR_end (external numbering )
	my($SEGLEN) = 8;
	my($PSGTHR) = 4.4;

    my($cr_leng, $poschg, $negchg, $ur_leng) = &regions;

	my $ur_leng2 = ($cr_leng+$ur_leng > 30) ? 30 - $cr_leng : $ur_leng;
	return -10000 if $cr_leng+$SEGLEN > $PsortSeq::len;

    my($ur_peak) = &hydseg2($SEGLEN, $cr_leng+1, $ur_leng2);


    if($OPT_V) {
		print "PSG:  a new signal peptide prediction method$endl";
	}
## ad hoc rules ##
	if($ur_leng >= 60 || ($poschg-$negchg) < 0) {
		$ur_peak = 0;
		print ">>> too long uncharged segment: $ur_leng$endl"
			if($ur_leng >= 60 && $OPT_V);
		printf ">>> negative net charge in the N-region: %3.0f$endl", 
			$poschg-$negchg if($poschg-$negchg < 0 && $OPT_V);
	}
    if($OPT_V) {
		printf("      N-region:  length %d;  pos.chg %d;  neg.chg %d$endl",
								$cr_leng, $poschg, $negchg);
		printf("      H-region:  length %d;  peak value %6.2f$endl",
              					$ur_leng, $ur_peak);
		printf("      PSG score: %6.2f$endl$endl", $ur_peak - $PSGTHR);
    }
	
	$ur_peak - $PSGTHR;
}

# return the postion of most C-term. charged residue
# (output 0 means 'no charge'),
# the charges within a region, and the ur_leng;
sub regions {
# calculate the CR length
    my($i, $aa);
    for( $i = $MaxPositionOf_CR_End; $i >= 0; --$i ) {
	$aa = $PsortSeq::seq[$i];
	last if($aa eq 'R' || $aa eq 'D' || $aa eq 'E' || $aa eq 'K');        }
    my $CRend = $i;

# count charged residues
    my($poschg, $negchg) = (0, 0);
    for ($i=0; $i<=$CRend; $i++) {
                $aa = $PsortSeq::seq[$i];
                if($aa eq 'R' || $aa eq 'K') {$poschg++}
                elsif($aa eq 'D' || $aa eq 'E') {$negchg++}
    }

# calculate the UR length
    my $ur_leng = 0;
    for( $i = $CRend+1; $i < $PsortSeq::len; $i++, $ur_leng++) {
                $aa = $PsortSeq::seq[$i];
                last if($aa eq 'R' || $aa eq 'D' || $aa eq 'E' || $aa eq 'K');        }
        #print substr($PsortSeq::seq, $start, $ur_leng), "  ";

        ($CRend+1, $poschg, $negchg, $ur_leng);
}

### most hydrophobic segment
### Args: SEGLENG, START, RANGE0
sub hydseg2 {
        my($segleng, $start, $ur_leng) = @_;
        die "*** too short sequence$endl" if $start + $segleng - 1 > $PsortSeq::len;
        $start--;  ### internal notation
        my($range) = ($ur_leng - $segleng < 0) ? 0 : ($ur_leng - $segleng);

        my(@aas) = split(//, substr($PsortSeq::seq, $start, $segleng+$range));
        my($i);
        my($sum) = 0;
        foreach $i (1..$segleng){
        }
        my($hmax, $imax) = ($sum, 0);
        foreach $i (1..$range) {
                $sum = $sum + $$ges{$aas[$i-1]} - $$ges{$aas[$i+$segleng-1]};
                if($sum > $hmax) {
                        $hmax = $sum;
                        $imax = $i;
                }
        }
        $hmax /= $segleng;
}


## von Heijne's method for signal seq prediction ##################
sub gvh {
    package GvH;

    my $organismType = $_[0];
    my ($i, $j, $end);

# parameter for prokaryotes
    my $prm_prk = [
    [ 1.14, 0.92, 0.92, 1.03, 0.63, 0.78, 0.45, 0.63, 0.78, 0.78,
      2.01,-0.47, 2.27, 1.73, 0.22 ],
    [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
     -3.58, 0.00,-3.58, 0.00, 0.00 ],
    [-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,-0.69,
     -3.58,-0.69,-3.58, 0.00, 1.39 ],
    [-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,-0.79,
     -3.58,-0.79,-3.58, 0.60, 1.29 ],
    [ 0.43, 1.12, 0.84, 1.12,-0.26,-0.26, 1.82,-0.26, 1.12,-0.26,
     -3.58, 1.68,-3.58,-0.26,-0.26 ],
    [ 0.39,-0.30,-0.30,-0.30, 0.11, 0.62,-0.30, 0.39,-0.30,-0.30,
     -3.58,-0.30,-0.30,-0.99,-0.99 ],
    [ 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22,
     -3.58, 2.17,-3.58, 0.22, 0.22 ],
    [ 0.57,-0.53, 1.08,-0.53, 1.08,-0.53,-0.53, 0.57,-0.53,-0.53,
     -3.58,-0.53,-3.58,-0.53, 0.16 ],
    [-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,-0.92,
     -3.58,-0.22,-3.58, 0.18,-0.92 ],
    [ 1.09, 1.40, 1.20, 1.09, 1.20, 1.57,-0.99,-0.99,-0.30,-0.30,
     -0.99,-0.30,-3.58,-0.99,-0.99 ],
    [ 0.51, 1.20, 0.51, 0.51, 1.61, 1.20, 1.61, 0.51, 0.51, 1.20,
     -3.58, 1.90,-3.58, 0.51, 0.51 ],
    [-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,-0.47,
     -3.58, 0.63,-3.58,-0.47, 0.92 ],
    [-0.53,-0.53,-0.53,-0.53,-0.53,-0.53, 0.16, 0.57, 1.08, 0.16,
     -3.58,-0.53,-3.58,-0.53, 1.08 ],
    [-0.34,-0.34,-0.34,-0.34,-0.34,-0.34,-0.34,-0.34, 0.36, 0.36,
     -3.58, 0.76,-3.58,-0.34,-0.34 ],
    [-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,-0.53,
     -3.58,-0.53,-3.58,-0.53,-0.53 ],
    [-0.96,-0.96,-0.96, 0.43, 0.43,-0.96, 0.65, 1.75, 0.65, 1.12,
      0.65,-0.26,-0.26,-0.96,-0.96 ],
    [-0.10,-0.79, 0.60,-0.10,-0.10,-0.10,-0.10,-0.10, 0.82,-0.79,
      0.31,-0.79,-0.79,-0.79,-0.10 ],
    [ 0.69, 1.03,-0.92, 0.18,-0.92, 0.47, 1.03,-0.92,-0.92, 0.47,
      0.18,-0.92,-3.58,-0.22,-0.92 ],
    [ 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92,
     -3.58, 0.92,-3.58, 0.92, 0.92 ],
    [-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26,-0.26, 0.84,
     -3.58,-0.26,-3.58,-0.26,-0.26 ]
   ];

# parameter for eukaryotes
    my $prm_euk = [
    [ 0.10,-0.11,-0.04, 0.03, 0.32, 0.22, 0.22, 0.16, 0.54, 0.03,
      1.18,-0.88, 1.71, 0.22,-0.88 ],
    [-0.41, 0.29, 0.69, 0.44, 0.69, 1.13, 0.29, 0.58, 0.11, 0.29,
      1.44,-0.41, 0.69, 0.58,-0.41 ],
    [-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-2.19,-0.58,-1.09,
     -5.08,-0.58,-5.08, 0.12, 0.21 ],
    [-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-2.30,-1.20,-0.36,
     -5.08,-0.36,-5.08, 0.26, 0.34 ],
    [ 0.84, 0.47, 0.68, 0.68, 0.07, 0.22, 1.17, 0.84,-0.34,-0.11,
     -5.08, 0.84,-5.08, 0.07,-0.34 ],
    [-1.11,-1.11,-1.39,-0.70,-1.39, 0.07,-1.39,-1.80, 0.45, 1.03,
     -0.88,-0.55, 1.17,-0.19,-0.55 ],
    [-1.22,-1.22,-1.22,-1.22,-1.22,-1.22,-1.22,-1.22, 0.39,-1.22,
     -5.08, 0.57,-5.08, 0.16,-0.53 ],
    [ 0.71, 0.71, 0.08,-0.21, 0.40,-0.39,-0.62, 0.08,-0.39,-2.00,
      0.30,-0.39,-5.08, 0.08,-0.06 ],
    [-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-2.42,-1.04,
     -5.08,-1.73,-5.08,-0.03,-0.23 ],
    [ 1.77, 1.73, 1.78, 1.88, 1.86, 1.31, 1.67, 1.40,-0.19, 0.64,
     -0.41, 0.50,-2.49,-0.41,-1.11 ],
    [-0.99, 0.11, 0.95, 0.39,-0.99, 0.80,-0.30,-0.30,-0.99,-0.99,
     -5.08,-0.99,-5.08,-0.99,-0.30 ],
    [-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-1.96,-0.86,-0.86,
     -5.08, 0.34,-5.08,-0.57,-0.01 ],
    [-1.31,-2.00,-1.31,-2.00,-2.00,-0.62,-2.00, 0.08, 0.99, 0.64,
     -5.08,-2.00,-0.90,-2.00, 1.09 ],
    [-1.84,-1.84,-1.84,-1.84,-1.84,-0.05,-1.84,-1.84, 0.46, 0.24,
     -5.08, 1.05,-0.74, 1.10, 0.46 ],
    [-1.34,-2.03,-2.03,-2.03,-2.03,-2.03,-2.03,-2.03,-0.08,-0.64,
     -5.08, 0.68,-5.08, 0.46, 0.17 ],
    [-0.24,-1.34,-0.35,-0.64, 0.13,-0.13, 0.27, 0.34, 0.82,-0.04,
      0.70, 0.40, 0.56, 0.27,-0.13 ],
    [-1.58, 0.03,-0.66,-0.89,-0.66, 0.29,-0.33,-0.33, 0.21,-0.48,
      0.56,-0.19,-0.48,-1.17, 0.03 ],
    [ 0.59, 0.81, 0.30, 0.48, 0.16, 0.30,-0.01, 0.89,-2.41, 0.08,
      1.06,-1.31,-5.08,-0.33, 0.43 ],
    [ 0.80, 0.51, 0.51,-0.59,-0.59, 0.11, 1.20, 0.51,-0.59, 0.51,
     -5.08, 1.61,-5.08, 0.11,-0.59 ],
    [-1.72,-1.72,-0.34,-1.72,-1.72,-1.72,-0.62,-1.72,-1.72,-1.03,
     -5.08,-0.11,-5.08,-1.72, 0.22 ]
   ];

## GvH numbering method
    my %aagvh = (
     'A'=> 0, 'C'=> 1, 'D'=> 2, 'E'=> 3, 'F'=> 4,
     'G'=> 5, 'H'=> 6, 'I'=> 7, 'K'=> 8, 'L'=> 9,
     'M'=>10, 'N'=>11, 'P'=>12, 'Q'=>13, 'R'=>14,
     'S'=>15, 'T'=>16, 'V'=>17, 'W'=>18, 'Y'=>19
    );

	my($imax, $smax) = (-1, -1000);
    	## position of max score  ## max score
    $end = ($PsortSeq::len - 15 > 48) ? 48 : ($PsortSeq::len - 15);
    my($score);

    foreach $i (0..$end) {
		$score = 0;
		foreach $j (0..14) {
	   		next unless defined($aagvh{$PsortSeq::seq[$i+$j]});
	    	if ($organismType eq "prokaryote") {
				$score += $prm_prk->[$aagvh{$PsortSeq::seq[$i+$j]}][$j];
	    	} elsif($organismType eq "eukaryote") {
				$score += $prm_euk->[$aagvh{$PsortSeq::seq[$i+$j]}][$j];
	    	} else {die "gvh: invalid argument"}
		}
#	printf "%2d: %5.2f$endl", $i+13, $score-7.5;
		if ($score > $smax) {
	    	$smax = $score;
	    	$imax = $i;
		}
    }
    my $gvhscore = $smax - 7.5;
    my $cleavsite = $imax + 13;


    if($OPT_V) {
	print "GvH:  von Heijne's method for signal seq. recognition$endl";
	printf "      GvH score (threshold: -2.1): %6.2f$endl", $gvhscore;
	printf "      possible cleavage site: between %2d and %2d$endl$endl", 
	$cleavsite, $cleavsite+1;
    }
    
    ($gvhscore, $cleavsite);
}


## Klein et al's ALOM ###############################################

# new version (2 threshold parameters)
sub alom2 {
    my($start, $threshold_loose, $threshold_strict) = @_;
	## threshold_loose (larger value) > threshold_strict (smaller)
	$start -= 1; ## internal residue numbering

    my ($A1, $A0) = (-9.02, 14.27);
    return (-1, 999, -1) if ($PsortSeq::len < 20 || $start > $PsortSeq::len-$LSEG);

    my($i, $j, $count, @hyd, $xmax, $sum, $imax, $hmax, $x0);
    undef %tms;
    ###  %tms is used by other routenes

    

    foreach $i ($start..$PsortSeq::len-1) {
	$hyd[$i] = $$kyteDoolittle{$PsortSeq::seq[$i]};
    }

	if($OPT_V) {
		print "ALOM: Klein et al's method for TM region allocation$endl";
		printf "      Init position for calculation: %d$endl", $start+1;
	}

    for($count=0; ; $count++) {

#### hydseg
		$sum = 0;
		$imax = $start;
		foreach $i ($start..$start+$LSEG-1) {
		    $sum += $hyd[$i];
		}
		$hmax = $sum;

		foreach $i ($start+1..$PsortSeq::len-$LSEG) {
			$sum = $sum - $hyd[$i-1] + $hyd[$i+$LSEG-1];
			if($sum > $hmax) {
				$hmax = $sum;
				$imax = $i;
			}
		}
		$hmax /= $LSEG;

#### allocation
		$x0 = $A1 * $hmax + $A0;
		$xmax = $x0 if $count==0;  ### maximaum hydrophobicity value

		last if $x0 > $threshold_loose;
		$tms{$imax+1} = $x0;

#### erase the segment
		foreach $i (0..$LSEG-1) {
	    	$hyd[$imax+$i] = -10;
		}
    }

#### check the strict threshold value
	printf "      Tentative number of TMS(s) for the threshold %4.1f:  %2d$endl", 
			$threshold_loose, $count if $OPT_V;
	if($count > 0) {
		my($count2) = 0;
		foreach (keys %tms) {
			$count2++ if $tms{$_} <= $threshold_strict;
		}
### caution!  this may ignore uncleaved signals
		if($count2 < 2) {
			foreach (keys %tms) {
				delete $tms{$_} if $tms{$_} > $threshold_strict;
			}
			$count = $count2;
			printf "      Number of TMS(s) for threshold %4.1f:  %2d$endl", 
				$threshold_loose, $count if $OPT_V;
		}
	} else {
		printf "      number of TMS(s) .. fixed$endl" if $OPT_V;
	}

#### output
	$mostn = -1;
	unless($count == 0) {
		foreach (sort {$a <=> $b} keys %tms) {
			if($OPT_V) {
	    		printf "      INTEGRAL    Likelihood =%6.2f", $tms{$_};
	    		printf "   Transmembrane %4d -%4d$endl", $_, $_+$LSEG-1;
			}
			$mostn = $_ if $mostn < 0;
		}
	}
	if($OPT_V) {
		printf("      PERIPHERAL  Likelihood =%6.2f (at %d)$endl", 
			$x0, $imax+1);
		printf "      ALOM score: %6.2f  (number of TMSs: %d)$endl$endl", 
			$xmax, $count;
	}

    ($count, $xmax, $mostn);
}

## Membrane Topology ####################################################

# Renamed some variables for clarity and fixed calculation bug. Paul Horton 2005/12/13.
sub mtop {
    my($signalSeqAnchorCentralPosition) = $_[0];
    my($start, $end);
    my $flankLen = 15;
    my $mtopScore;
    if( $PsortSeq::len < 20 
	|| $signalSeqAnchorCentralPosition < 1 
	|| $signalSeqAnchorCentralPosition > $PsortSeq::len-1){
	return -100;
    }

    # set $start to first charged residue before $signalSeqAnchorCentralPosition
    for ($start = $signalSeqAnchorCentralPosition; $start >= 0; $start--) {
	my $curResidue = $PsortSeq::seq[$start];
	last if ( $curResidue eq 'R' || $curResidue eq 'D' || $curResidue eq 'E'
		  || $curResidue eq 'H' || $curResidue eq 'K');
    }
    
    # set $end to first charged residue after $signalSeqAnchorCentralPosition
    for($end = $signalSeqAnchorCentralPosition; $end < $PsortSeq::len; $end++){
	my $curResidue = $PsortSeq::seq[$end];
	last if ( $curResidue eq 'R' || $curResidue eq 'D' || $curResidue eq 'E'
		  || $curResidue eq 'H' || $curResidue eq 'K');
    }

    # charge difference
    my $NterminalFlankCharge = 0.0;
    for( my $curPos = $start; $curPos >= 0 && $start - $curPos < $flankLen; --$curPos ){
	my $curResidue = $PsortSeq::seq[$curPos];
	if($curResidue eq 'D' || $curResidue eq 'E') {
	    $NterminalFlankCharge -= 1.0;
	} elsif($curResidue eq 'R' || $curResidue eq 'K') {
	    $NterminalFlankCharge += 1.0;
	} elsif($curResidue eq 'H') {
	    $NterminalFlankCharge += 0.5;
	}
    }
    if( $start <= $flankLen ){
	# if N-terminus included in N-terminal side flank, add 1.0 for the N-terminal amino group.
	$NterminalFlankCharge += 1.0;
    }

    my $CterminalFlankCharge = 0.0;
    for( my $curPos = $end; $curPos < $PsortSeq::len && $curPos - $end < $flankLen; ++$curPos ){
	my $curResidue = $PsortSeq::seq[$curPos];
	if($curResidue eq 'D' || $curResidue eq 'E') {
	    $CterminalFlankCharge -= 1.0;
	} elsif($curResidue eq 'R' || $curResidue eq 'K') {
	    $CterminalFlankCharge += 1.0;
	} elsif($curResidue eq 'H') {
	    $CterminalFlankCharge += 0.5;
	}
    }

    
    $mtopScore = $CterminalFlankCharge - $NterminalFlankCharge;
    if($OPT_V) {
	print "MTOP: Prediction of membrane topology (Hartmann et al.)$endl";
	print ' 'x6, "Center position for calculation: $signalSeqAnchorCentralPosition$endl";
	printf "      Charge difference: %4.1f   C(%4.1f) - N(%4.1f)$endl", 
	$mtopScore, $CterminalFlankCharge, $NterminalFlankCharge;
	print ' ' x 6;
	($mtopScore>0) ? print "C > N: C-terminal side will be inside$endl$endl" :
	    print "N >= C: N-terminal side will be inside$endl$endl";
    }
    return $mtopScore;
}

## assignment of membrane topology

##   output: mtype, start, end
##
###  caution: uses $mostn
sub mtype_assign {
    my($tmsnum, $sig, $mtopscr) = @_;

    my($pos, $c);

    if($tmsnum == 0) {
	('__', 0, 0);
    } elsif($tmsnum == 1) {
	if($sig eq 'cleavable') {
	    $pos = (keys %tms)[0] + $LSEG;
	    $pos=$PsortSeq::len if $pos > $PsortSeq::len;  ##### ad hoc
	    ('1a', $pos, $PsortSeq::len);
	} elsif(($mostn / $PsortSeq::len) > 0.8) {  ### exceptional topology
	    print ">>> Single TMS is located near the C-terminus$endl$endl"
				if $OPT_V;
	    $pos = (keys %tms)[0] - 1;
	    ('Nt', 1, $pos);
	} else {
	    if($mtopscr > 0) {
		('1b', $mostn, $PsortSeq::len);
	    } else {
		('2 ', 1, $mostn);
	    }
	}
    } else {
	if($mtopscr > 0) {
	    ('3b', 0, 0);
	} else {
	    ('3a', 0, 0);
	}
    }
}


1;
