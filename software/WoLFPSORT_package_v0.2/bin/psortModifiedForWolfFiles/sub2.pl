# Author: Kenta Nakai
# Maintainer: Paul Horton (in the context of WoLF PSORT)
# Date of last CVS update: $Date: 2006/08/27 11:02:13 $
#
############################################################
##  PSORT II:  subroutine package 2
##                        copyright: Kenta Nakai 1996, 1997
##                               last update: May 3, 1997
############################################################

# 2006/07/12 added "our" declaration of global variables.

# This version is for PSORT+WoLF, slightly modified from original May 2004.
# The variable $endl was added by Paul Horton to allow the end of line string
# to be changed (e.g. to "<br>\n").
BEGIN{  
    use FindBin;
    push @INC, "$FindBin::Bin/psortModules";
}
use PsortAminoAcidIndexer;
use PsortAminoAcidIndex;
#our ( $OPT_V, $endl, $seq, $seqLen );
our ( $OPT_V, $endl );
use strict;
package sub2;

## Sequence Numbering convention:
##   argument: 1, 2, 3,...   internal: 0, 1, 2,...

my $LSEG2 = 8;


## mitochondrial #############################################

sub mitdisc {
	my( $m75, $m95, $limit, $score);
	my @freq;

	foreach my $i (0..PsortAminoAcidIndexer::size()-1) {$freq[$i] = 0}
	my $range = &getrange($PsortSeq::seq);
	foreach my $j (0..$range-1) {
	    my $index = PsortAminoAcidIndexer::toIndex( $PsortSeq::seq[$j] );
	    $freq[$index]++;
	}
	$limit = ($range>$LSEG2) ? $LSEG2 : $range;
	my $sequenceSegment = substr( $PsortSeq::seq, 0, $limit );
	my @sequenceSegment = split( //, $sequenceSegment );
	$m95 = &hmom(95, \@sequenceSegment );
        $m75 = &hmom(75, \@sequenceSegment );

	my(@a) = 
		( 2.0102 - 1.0973,  ## 0: R
		  0.3922 - 0.2538,  ## 1: M75
		  0.4737 - 0.3566,  ## 2: M95
		 -0.7417 + 0.2655,  ## 3: G 
		  9.2710 -11.3204,  ## 4: AC 
		  0.7522 - 0.4503,  ## 5: ST 
		-17.5993 +13.0901 ); ##6: const

	$score = $a[0] * $freq[1]
		   + $a[1] * $m75
		   + $a[2] * $m95
		   + $a[3] * $freq[7]
		   + $a[4] * ($freq[3]+$freq[6])
		   + $a[5] * ($freq[15]+$freq[16])
		   + $a[6];

	if($OPT_V) {
		print "MITDISC: discrimination of mitochondrial targeting seq$endl";
		printf "      R content:   %5d       Hyd Moment(75): %5.2f$endl", 
						$freq[1], $m75;
		printf "      Hyd Moment(95): %5.2f    G content:   %5d$endl", 
						$m95, $freq[7];
		printf "      D/E content: %5d       S/T content: %5d$endl",
						$freq[3]+$freq[6], $freq[15]+$freq[16];
		printf "      Score: %5.2f$endl$endl", $score;
	}
	$score;
}

sub hmom {
    my ($deg, $sequenceSegmentRef) = @_;
    my $scr;
    my($radianAngle) = 3.14159 * $deg / 180.0;
    my ($sumSin, $sumCos) = (0, 0);
    my $sequenceSegmentLength = scalar @$sequenceSegmentRef;

    my $ges = PsortAminoAcidIndex->new( "Engelman GES 1986" );

    foreach my $i (0..$sequenceSegmentLength-1) {
	$sumSin += $$ges{ $$sequenceSegmentRef[$i] } * sin( $radianAngle * $i );
	$sumCos += $$ges{ $$sequenceSegmentRef[$i] } * cos( $radianAngle * $i );
    }
    $scr = sqrt($sumSin ** 2 + $sumCos ** 2) / $LSEG2;
    # maybe this should be dividing by sequenceSegmentLength instead? P.H. 2006/4/28.
}

sub getrange {
        my($sq) = $_[0];
#		print $sq, "$endl";

        my($POS1) = 15;
        my($pos, $i, $c); ## probable transit peptide 1-$pos
        my($flg) = "not found"; ## acidic residue before $POS1
        my($len) = length($sq);
        if($len < $POS1) {return 0}
        for($i=1; $i<$POS1; $i++) {
                $c = substr($sq, $i, 1);
                if($c eq 'D' || $c eq 'E') {
                        $flg = "found";
                        $pos = $i+1; # from the next residue
                        last;
                }
        }
        if($flg eq "not found") {$pos = $POS1}
        for($i=$pos; $i<$len; $i++) {
                $c = substr($sq, $i, 1);
                if($c eq 'D' || $c eq 'E') {
                        return $i+1; ## pos is counted from 1
                }
        }
        $len;  # if no acidic resudue at all
}

## gavel et al's method
sub gavel {
	my($ipos, $isite, $motif);
	print "Gavel: prediction of cleavage sites for mitochondrial preseq$endl"
			if $OPT_V;
	if(($ipos = &r3) > 0) {
		$isite = $ipos + 4;
		$motif = 'R-3';
		if($OPT_V) {
			print "      R-3 motif at $isite  ";
			print substr($PsortSeq::seq, $ipos-1, 4),"|",
				substr($PsortSeq::seq, $ipos+3, 1),"$endl$endl";
		}
	} elsif (($ipos = &r10) > 0) {
		$isite = $ipos + 11;
		$motif = 'R-10';
		if($OPT_V) {
			print "      R-10 motif at $isite  ";
			print substr($PsortSeq::seq, $ipos-1, 3)," ",
				substr($PsortSeq::seq, $ipos+2, 2),"$endl$endl";
		}
	} elsif (($ipos = &r2) > 0) {
		$motif = 'R-2';
		$isite = $ipos + 11;
		if($OPT_V) {
			print "      R-2 motif at $isite  ";
			print substr($PsortSeq::seq, $ipos-1, 3),"|",
					substr($PsortSeq::seq, $ipos+2, 2),"$endl$endl";
		}
	} else {
		$isite = 0;
		$motif = '___';
		if($OPT_V) {
			print "      cleavage site motif not found$endl$endl";
		}
	}
	($isite, $motif);
}

sub abs {
	($_[0]>=0) ? $_[0] : (-1 * $_[0]);
}

## RXY(S/A) motif
sub r3 {
	my($max) = ($PsortSeq::len < 67) ? ($PsortSeq::len - 10) : 67;
	my($min) = 11;  ## 12-1
	my($ipos) = -1000;
	foreach my $i ($min..$max-3) {  ## $max-1 -> $max-3
		if($PsortSeq::seq[$i] eq 'R' && $PsortSeq::seq[$i+2] eq 'Y'
			&& ($PsortSeq::seq[$i+3] eq 'S' || $PsortSeq::seq[$i+3] eq 'A')) {
			$ipos = $i if &abs($i-23) < &abs($ipos-23);
		}
	}
	$ipos;
}

## RXFS motif
sub r10 {
	my($max) = ($PsortSeq::len < 60) ? ($PsortSeq::len-10) : 60;
	my($min) = 4;  ## 5-1
	my($ipos) = -1000;
	foreach my $i ($min..$max-3) {  ## $max-1 -> $max-3
		if($PsortSeq::seq[$i] eq 'R' && $PsortSeq::seq[$i+2] eq 'F' && $PsortSeq::seq[$i+3] eq 'S') {
			$ipos = $i if &abs($i-23) < &abs($ipos-23);
		}
	}
	$ipos;
}

## transition point in terms of negative charges
sub r2 {
	my($i, $j, $ipos);
	my($k) = -1000;
	LOOP: foreach $i (0..$PsortSeq::len-13) {
		if($PsortSeq::seq[$i] eq 'D' || $PsortSeq::seq[$i] eq 'E') {
			foreach $j ($i+1..$i+12) {
				if($PsortSeq::seq[$j] eq 'D' || $PsortSeq::seq[$j] eq 'E') {
					$k = 0;  ### found
					$ipos = $i;
					last LOOP;
				}
			}
		}
	}
	if($k == 0) {
		for ($k = $ipos-2; $k > 0; $k--) {
			last if $PsortSeq::seq[$k] eq 'R';
		}
	}
	return $k;
}

# nuclear ##############################################################

sub nucdisc {
	my($i, $j, $score);

	print "NUCDISC: discrimination of nuclear localization signals$endl"
				if $OPT_V;

### new parameter 
#$start = (times)[0];
	$score = 0;
	$score += ( 0.0901 - 0.0274) * &nls1;
	$score += ( 0.1648 - 0.0786) * &nls2;
	$score += ( 0.1063 - 0.0241) * &bipartite;
	$score += ( 1.1162 - 0.1665) * &nucaa;
	$score += (-1.2642 + 0.7904);
#$end = (times)[0];
#printf "(old)%.2f\t", $end - $start;

	printf "      NLS Score: %5.2f$endl$endl", $score if $OPT_V;
	$score;
}

# nls1: 4 residue pattern
sub nls1 {
	my($i, $j, $scr, $c, $iscr);
	my($nbas, $np, $nh, $flg) = (0, 0, 0, 'none');
# init value
	foreach $j (0..3) {
		$c = $PsortSeq::seq[$j];
		if($c eq 'R' || $c eq 'K') {$nbas++;}
		elsif($c eq 'P') {$np++;}
		elsif($c eq 'H') {$nh++;}
	}
	if($nbas == 4) {$scr = 5}
	elsif($nbas == 3 && $np == 1) {$scr = 4}
	elsif($nbas == 3 && $nh == 1) {$scr = 3}
	else {$scr =  0}
	if($scr > 0 && $OPT_V) {
		$flg = 'found';
		printf "      pat1: %s at    1$endl", substr($PsortSeq::seq, 0, 4);
	}

	foreach $i (0..$PsortSeq::len-5) {
		$c = $PsortSeq::seq[$i+4];
		if($c eq 'R' || $c eq 'K') {$nbas++;}
		elsif($c eq 'P') {$np++;}
		elsif($c eq 'H') {$nh++;}
		$c = $PsortSeq::seq[$i];
		if($c eq 'R' || $c eq 'K') {$nbas--;}
		elsif($c eq 'P') {$np--;}
		elsif($c eq 'H') {$nh--;}

		if($nbas == 4) {$iscr = 5}
		elsif($nbas == 3 && $np == 1) {$iscr = 4}
		elsif($nbas == 3 && $nh == 1) {$iscr = 3}
		else {$iscr = 0}
		if($iscr > 0) {
			$scr += $iscr;
			$flg = 'found' if $flg eq 'none';
			printf "      pat4: %s (%1d) at %4d$endl", 
				substr($PsortSeq::seq, $i+1, 4), $iscr, $i+2 if $OPT_V;
		}
	}
	printf "      pat4: none$endl" if $OPT_V && ($flg eq 'none');
	
	$scr;
}

# nls2: (max) 7 residue pattern
sub nls2 {
	my($i, $iscr);
	my($scr, $flg) = (0, 'none');
	foreach $i (0..$PsortSeq::len-8) {
		if($PsortSeq::seq[$i] eq 'P') {
			$iscr = &scr7($i);
			if($iscr != -1) {
				$scr += $iscr;
				$flg = 'found' if $flg eq 'none';
				printf "      pat7: %s (%d) at %4d$endl", 
					substr($PsortSeq::seq, $i, 7), $iscr, $i+1 if $OPT_V;
			}
		}
	}
	printf "      pat7: none$endl" if($OPT_V && ($flg eq 'none'));
	$scr;
}

sub scr7 {
	my($i) = $_[0];
	return -1 unless $PsortSeq::seq[$i] eq 'P';
	my($max) = -1;
	my($j, $k, $nbas);
	foreach $k (0..2) {
		$nbas = 0;
		foreach $j (1..4) {
			$nbas++ if $PsortSeq::seq[$i+$j+$k]eq'R'||$PsortSeq::seq[$i+$j+$k]eq'K'; 
		}
		if($nbas == 4) {$max = 5}
		elsif($nbas == 3) {$max = ($max > 5-$k) ? $max : (5-$k)}
	}
	$max;
}

# robbins: bipartite NLS (Robbins ,, Dingwall)

sub bipartite {
    my $cnt;
    my $scr = 0;
    my $numBipartiteFound = 0;
    foreach my $i (0..$PsortSeq::len-17) {
	next unless (substr($PsortSeq::seq, $i, 2) =~ /[RK][RK]/);
	$cnt = 0;
	foreach my $j ($i+12..$i+16) {
	    $cnt++ if $PsortSeq::seq[$j] eq 'R' || $PsortSeq::seq[$j] eq 'K';
	}
# *************************** LOOK HERE!!! *****************************
# adding cnt to scr twice looks like it is almost certainly a bug.
# However changing it at this instant would require changing the dataset
# which I do not want to do. Remember to get rid of this in the future!
# P.H. 2006/7/3.
	$scr += $cnt if $cnt >= 3; # bug!!!
# # ********************************************************************
	if($cnt >= 3) {
	    $scr += $cnt;
	    ++$numBipartiteFound;
	    if( $OPT_V ){
		printf "      bipartite: %s at %4d$endl", substr($PsortSeq::seq, $i, 17), $i+1;
	    }
	}
    }
    if( ($OPT_V) && (!$numBipartiteFound) ){
	printf "      bipartite: none$endl";
    }
    $scr;
}

# nucaac
sub nucaa {
	my($cnt) = 0;
	my($i);
	foreach $i (0..$PsortSeq::len-1) {
		if($PsortSeq::seq[$i]eq'R' || $PsortSeq::seq[$i]eq'K') {
			$cnt++;
		}
	}
	$cnt /= $PsortSeq::len;
	printf "      content of basic residues: %5.1f%%$endl", $cnt * 100
			if $OPT_V;
	if($cnt > 0.2) {
		int($cnt * 10 - 1);
	} else {
		0;
	}
}

# RNA-binding motif #######################################################

sub rnp1 {
	my(@match);
	print "RNA-binding motif:" if $OPT_V;
	@match = ($PsortSeq::seq =~ /[RK]G[^EDRKHPCG][AGSCI][FY][LIVA].[FYM]/g);
	if($OPT_V) {
		&report_motif(0, @match);
		print "$endl";
	}
	($#match < 0) ? 0 : ($#match+1);
}

# actinin-type actin-binding ##############################################

sub actin {
	my(@match);
	my($hit) = 0;
	if($OPT_V) {
		print "Actinin-type actin-binding motif:$endl";
		print "      type 1:";
	}
	@match = ($PsortSeq::seq =~ /[EQ]..[ATV]F..W.N/g);
	if($OPT_V) {
		&report_motif(0, @match);
	}
	$hit += scalar(@match);
	print "      type 2:" if $OPT_V;
	@match = ($PsortSeq::seq =~ /[LIVM].[SGN][LIVM][DAGHE][SAG].[DEAG][LIVM].[DEAG]....[LIVM].L[SAG][LIVM][LIVM]W.[LIVM][LIVM]/g);
	if($OPT_V) {
		&report_motif(0, @match);
	}
	print "$endl" if $OPT_V;
	$hit += scalar(@match);
}

# lumen pf ER ##############################################################

sub hdel {
	my($seg) = substr($PsortSeq::seq, $PsortSeq::len-4);
	print "KDEL: ER retention motif in the C-terminus: " if $OPT_V;

	if ($seg =~ /HDEL$/ || $seg =~ /KDEL$/) {
		print "$seg$endl$endl" if $OPT_V;
		return 1;
	} else {
		print "none$endl$endl" if $OPT_V;
		return 0;
	}
}

# peroxisome ##############################################################

sub pts1 {
    my($scr, $pat);
    print "SKL: peroxisomal targeting signal in the C-terminus: " if $OPT_V;
		
    if($PsortSeq::seq    =~ /SKL$/) {$pat = 'SKL'; $scr = 10/12}
    elsif($PsortSeq::seq =~ /SKF$/) {$pat = 'SKF'; $scr = 1/2}
    elsif($PsortSeq::seq =~ /AKL$/) {$pat = 'AKL'; $scr = 1/2}
    elsif($PsortSeq::seq =~ /([SAGCN][RKH][LIVMAF])$/) {$pat = $1; $scr = 1/4}
    else {$scr = 0}
    if($OPT_V) {
	if($scr > 0) {
	    print "$pat$endl$endl";
	} else {
	    print "none$endl$endl";
	}
    }
    $scr;
}

sub pts2 {
	my($scr, @match);
	print "PTS2: 2nd peroxisomal targeting signal: "
		if $OPT_V;
	@match = ($PsortSeq::seq =~ /[RK][LI].....[HQ]L/g);
	if($OPT_V) {
		&report_motif(0, @match);
		print "$endl";
	}
	($#match < 0) ? 0 : ($#match+1);
}

# vacuolar ##############################################################

sub vaccalc {
	my(@match);
	print "VAC: possible vacuolar targeting motif:" if $OPT_V;
	@match = ($PsortSeq::seq =~ /[TIK]LP[NKI]/g);
	if($OPT_V) {
		&report_motif(0, @match);
		print "$endl";
	}
	($#match < 0) ? 0 : ($#match+1);
}

# N-myristoyl ############################################################

sub nmyr {
	my($pat, $i);
	my($scr, $num_k) = (0, 0);
	print "NMYR: N-myristoylation pattern : " if $OPT_V;
	if($PsortSeq::seq =~ /^(M?G[^EDRKHPFYW]..[STAGCN][^P])/) {
		$pat = $1;
		$scr++;
		print $pat, "$endl" if $OPT_V;
		if($PsortSeq::seq =~ /^M?GC/) {
			$scr++;
			print "      3rd aa is cysteine (may be palmitylated)$endl"
				if $OPT_V;
		}
		foreach $i (3..9) {
			$num_k++ if $PsortSeq::seq[$i] eq 'K'
		}
		if($num_k >= 2) {
			$scr++;
			print "      additional alternating lysine motif$endl"
				if $OPT_V;
		}
		print "$endl" if $OPT_V;
	} else {
		print "none$endl$endl" if $OPT_V;
	}
	$scr;
}

# isoprenyl ##############################################################

sub isoprenyl {
	print "Farnesylation/Geranylgeranylation motif: " if $OPT_V;
	if($PsortSeq::seq =~ /(C[^DENQ][LIVM].)$/) {
		if($OPT_V) {
			print "      CaaX motif in the C-terminus: $1$endl";
			print "         if X is S, A, or M, it will be farnesylated$endl";
			print "         otherwise, it will be geranylgeranylated$endl$endl";
		}
		2;
	} elsif($PsortSeq::seq =~ /(C.C)$/) {
		print "      CXC motif in the C-terminus: $1$endl$endl" 
			if $OPT_V;
		1;
	} elsif($PsortSeq::seq =~ /(CC..)$/) {
		print "      CC motif near the C-terminus: $1$endl$endl" 
			if $OPT_V;
		1;
	} else {
		print "none$endl$endl" if $OPT_V;
		0;
	}
}

# ER membrane ##############################################################

sub erm {
    my($mtype) = $_[0];
    my($count, $scr) = (0, 0);
    my(@pat);
    print "ER Membrane Retention Signals: " if $OPT_V;
    my($nseg, $cseg) = (substr($PsortSeq::seq, 1, 4), substr($PsortSeq::seq, $PsortSeq::len-5, 4));

    @pat = ($nseg =~ /R/g);
    $count = scalar(@pat);
    if($count != 0) {
	print "$endl      XXRR-like motif in the N-terminus: $nseg$endl$endl"
	    if $OPT_V;
	$scr = $count;
	$scr += 2 if $mtype eq '2 ' || $mtype eq 'Nt';
    }

    @pat = ($cseg =~ /K/g);
    $count = scalar(@pat);
    if($count != 0) {
	print "$endl      KKXX-like motif in the C-terminus: $cseg$endl$endl"
	    if $OPT_V;
	$scr += $count;
	$scr += 2 if $mtype eq '1a' || $mtype eq '1b';
    } else {
	print "none$endl$endl" if $OPT_V;
    }
    $scr;
}

# YQRL motif ##############################################################

sub yqrl {
    my($tmsnum, $start) = @_;

    my(@match);
    my($scr) = 0;
    print "memYQRL: transport motif from cell surface to Golgi:" if $OPT_V;
    if($tmsnum == 0) {
	print " none$endl$endl" if $OPT_V;
    } elsif($tmsnum == 1) {
#	@match = (substr($PsortSeq::seq, $start-1, $end-$start+1) =~ /YQRL/g);
	@match = ($PsortSeq::seq =~ /YQRL/g);
	$scr += scalar(@match) + 3;
	if($OPT_V) {
	    &report_motif($start-1, @match);
	    print "$endl";
	}
    } else {
	@match = ($PsortSeq::seq =~ /YQRL/g);
	$scr += scalar(@match);
	if($OPT_V) {
	    &report_motif($start-1, @match);
	    print "$endl";
	}
    }
    $scr;
}

# tyrosine(s) in the tail #################################################

sub tyros {
    my($tmsnum, $start, $end) = @_;
    my(@ylist) = ();
    my($i);
    print "Tyrosines in the tail: " if $OPT_V;
    if($tmsnum != 1) {
	print " none$endl$endl" if $OPT_V;
	return 0;
    } elsif($end - $start + 1 > 50) {
	print " too long tail$endl$endl" if $OPT_V;
	return 0;
    } else {
	for($i=$start-1; $i<=$end-1; $i++) {
	    push(@ylist, $i) if $PsortSeq::seq[$i] eq 'Y';
	}
	if(scalar(@ylist) == 0) {
	    print " none$endl$endl" if $OPT_V;
	    return 0;
	}
	print join(",", @ylist), "$endl$endl" if $OPT_V;
	return 10 * scalar(@ylist) / ($end - $start + 1);
    }
}

# dileucine motif ########################################################

sub dileu {
    my($tmsnum, $start, $end) = @_;
    my(@match);
    my($scr) = 0;
    print "Dileucine motif in the tail:" if $OPT_V;
    if($tmsnum != 1) {
	print " none$endl$endl" if $OPT_V;
    } else {
	@match = (substr($PsortSeq::seq, $start-1, $end-$start+1) =~ /LL/g);
	$scr = 10 * scalar(@match) / ($end - $start + 1);
	if($OPT_V) {
	    &report_motif($start-1, @match);
	    print "$endl";
	}
    }
    $scr;
}

##########################################################################
## must be called when $OPT_V
sub report_motif {
	my($init, @pat) = @_;
##	print "### ", join(" ", @pat), "$endl";
	if($#pat < 0) {
		print " none$endl";
		return;
	}
	my($pos);
	print " found$endl";
	foreach (@pat) {
		$pos = index($PsortSeq::seq, $_, $init);
		$init = $pos + 1;
		printf "      %s at %d$endl", $_, $pos+1;
	}
}
	
1;
