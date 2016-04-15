#!/usr/bin/perl -w
#  Author: Paul B. Horton
#  Organization: Computational Biology Research Center, AIST, Japan
#  Creation Date: 2005.10.5
#  Last Modification Date: $Date: 2006/08/27 10:59:47 $
#  Copyright: Paul Horton
#  
#  Description: Read fasta format streams one record at a time.
#
#  History: This file branched off from the fasta reader stream
#           made by Paul Horton for the planned fastafmt package.
#
use Carp;
use Getopt::Long qw(:config posix_default bundling permute pass_through);
use strict;

package PsortFastaReader;
$PsortFastaReader::VERSION = 1.00;

# Maps user option names to their default values.
# Also used by _u() to test for the existence of an option name.
my %userOptionDefault = (
			 'add-asterisk' => 0,
			 'canonicalize' => 0,
			 'comment-regex' => "^\#",
			 'input-delimiter' => '\s+',
			 'no-comment' => 0,
			 'silent' => 0,
			 'strip-asterisk' => 0,
			 'strip-comments' => 0,
			 'strip-white-space-before-id' => 0,
			 'tab-delimited' => 0,
			 'warn' => 0,
			 'warn-level' => 0
 );

# user exposed options names. Passes though known keys.
# Croak on unknown names. Used to detect typing errors.
# Mnemonic "u" for "user".
sub _u{
    my $optionName = shift;
    defined( $userOptionDefault{$optionName} ) && return $optionName;
    Carp::croak( "Catastrophic Error: program tried to access unknown userOption \"$optionName\".
Might need to send a bug report" );
}

# Maps data variables to their initial values
# Also Used by _l() to test for the existence of an option name.
my %localVariableInitValue = (
			       _streamHeaderLines      => [],
			       _prevRecordLine         => "",
			       _atEndOfStream          => 0,
			       _recordLines            => [],
			       _recordWarnings         => [] );


# user exposed options names. Passes though known keys.
# Croak on unknown names.
# Mnemonic "l" for "local"
sub _l{
    my $variableName = shift;
    defined( $localVariableInitValue{$variableName} ) && return $variableName;
    Carp::croak( "Catastrophic Error: program tried to access unknown userOption \"$variableName\".
Might need to send a bug report" );
}

my %userOptionToShortCut =
    ( _u("warn")         => "w",
      _u("tab-delimited") => "|t" );

# If appropriate, add short-cut to option name to return
# the fullName|shortcutName needed by getOptions. The
# specifier for flags taking arguments, e.g. "=s", are not added
sub _toSpec{
    my $bareOpt = shift;
    my $optSpec;
    _u($bareOpt); # check for typing errors, etc.
    if(  defined( $userOptionToShortCut{$bareOpt} )  ){
	$optSpec = $bareOpt . "|" . $userOptionToShortCut{$bareOpt};
    }else{
	$optSpec = $bareOpt;
    }
    return $optSpec;
}


# Explanations should be *very* terse. Try to fit all on one line with
# an 80 character terminal.
my %optionOneLiner = 
    ( _toSpec("add-asterisk" )     => "make each printed record body end in \"*\"",
      _toSpec( "canonicalize" )    => "output records in canonical form",
      _toSpec( "comment-regex" )   => "regex used to define comment lines",
      _toSpec( "input-delimiter" ) => "character for delimiting fields in head lines",
      _toSpec( "no-comment" )      => "never consider any line a comment line",
      _toSpec( "silent" )          => "do not print warnings",
      _toSpec( "strip-asterisk" )  => "at end of record, do not print a trailing \"*\"",
      _toSpec( "strip-comments" )  => "do not output comment lines",
      _toSpec( "strip-white-space-before-id" ) => "remove any white space after the record head \">\"",
      _toSpec( "warn" )            => "be verbose about warnings",
      _toSpec( "tab-delimited" )   => "Only treat tabs as white space"
      );


# croak if %userOptionDefault contains options not defined in %optionOneLiner.
sub _checkOptionOneLinerCompleteness{
    
    for my $userOption (keys %userOptionDefault){
	next if( $userOption =~ /^_/ );  # options starting with "_" should not user exposed documentation.
	my $optionWithShortCut = 
	    defined($userOptionToShortCut{$userOption}) ? $userOptionToShortCut{$userOption} : $userOption;
	unless(  defined( $optionOneLiner{$optionWithShortCut} )  ){
	    Carp::croak( "Catastrophic Error: one liner documentation not defined for option $optionWithShortCut" );
	  }
    }
}

sub getOneLiners{
    return %optionOneLiner; # stopped programming here...
}

sub new{
    my( $class, %option ) = @_;
    for my $key (keys %option){
	defined $userOptionDefault{$key} 
	|| defined $localVariableInitValue{$key}
	|| Carp::croak( "unknown option: \"$key\"" );
    }
    my %data = ( %localVariableInitValue, %userOptionDefault, %option );
    my $dataRef = \%data;
    bless $dataRef, $class;
    return $dataRef;
}


sub getOptions{
    my $this = shift;

    defined( $this->{_recordWarnings} ) || die;
    my %spec =
  (
   _toSpec( "add-asterisk" ) =>
   sub{ $this->_set( "add-asterisk", 1 ) },
	 
   _toSpec( "canonicalize" ) =>
   sub{ $this->_set( "canonicalize", 1 ) },

   _toSpec( "comment-regex" )."=s" =>
   sub{ $this->_set( "comment-regex", $_[1] ) },

   _toSpec( "input-delimiter", )."=s" =>
   sub{ $this->_set( "inputDelimiter", $_[1] ) },

   _toSpec( "no-comment" ) =>
   sub{ $this->_set( "no-comment", 1 ) },

   _toSpec( "silent" ) =>
   sub{ $this->_set( "warn-level",   -1 ) },

   _toSpec( "strip-asterisk" ) =>
   sub{ $this->_set( "strip-asterisk", 1 ) },

   _toSpec( "strip-comments" ) =>
   sub{ $this->_set( "strip-comments", 0 ) },

   _toSpec( "strip-white-space-before-id" ) =>
   sub{ $this->_set( "strip-white-space-before-id", 1 ) },

   _toSpec( "tab-delimited" ) =>
   sub{ $this->_set( "tab-delimited", 1 );
	$this->_set( "input-delimiter", "\t+" ) },

   _toSpec( "warn" ) =>
   sub{ $this->_set( "warn-level", 1 ) }

); # end %spec initialization.

   my $getOptionsRetval = ::GetOptions( %spec );

   $this->_processOptions();
}

sub _processOptions(){
    my $this = shift;
    my $conflictingOptionsErrorMessage 
	= "Error: conflicting options \"%s\" and \"%s\"\n";
    if(  ( $this->get("add-asterisk") ) && ( $this->get("strip-asterisk") )  ){
	die(  sprintf( $conflictingOptionsErrorMessage, _u("add-asterisk"), u_("strip-asterisk") )  );
    }
    if( $this->get("canonicalize") ){
	if( $this->get("strip-asterisk") ){
	    die(  sprintf( $conflictingOptionsErrorMessage, _u("canonicalize"), u_("strip-asterisk") )  )
	}
	$this->_set( "add-asterisk", 1 );
	$this->_set( "strip-white-space-before-id", 1 );
	$this->_set( "strip-comments", 1 );
    }
}

sub optionNames{
    my $this = shift;
    return grep /^[^_]/, keys %{$this};
}

sub printOptionValues{
    my $this = shift;
    my $ostream = shift || *STDOUT;

    print "###value of add-asterisk is: ", $this->get("add-asterisk"), "\n";
    _prettyPrintHashValues( $ostream, $this, $this->optionNames() );
}


##################################################
# Begin main methods
##################################################
sub printRecordWarnings{
    my $this = shift;
    my $ostream = shift || *STDERR;

    my $warningRef = $this->get("_recordWarnings");
    for my $warning (@$warningRef){
	print $ostream "$warning\n";
    }
}
	
sub recordId{
    my $this = shift;

    return $this->recordField( 0 );
}

sub recordHead{
    my $this = shift;
    my $recordLinesRef = $this->get("_recordLines");
    return ""  unless( @$recordLinesRef );
    my $head = $$recordLinesRef[0];
    chomp $head;
    return $head;
}

sub recordDescription{ # head without the id
    my $this = shift;
    my $recordLinesRef = $this->get("_recordLines");
    return ""  unless( @$recordLinesRef );
    my $head = $$recordLinesRef[0];
    chomp $head;
    my $del = $this->get( 'input-delimiter' );
    if( $head =~ /$del(.*)$/ ){
	return $1;
    }
    return "";
}

sub recordField{
    my( $this, $fieldNum ) = @_;
    my $head = $this->recordHead();
    $head =~ s/^>\s*//;
    my $del = $this->get( 'input-delimiter' );
    my @headFields = split( /$del/, $head );
    return $headFields[$fieldNum];
}

sub recordBody{
    my $this = shift;
    my $recordLinesRef = $this->{_recordLines};
    return @$recordLinesRef[1 .. $#$recordLinesRef];
}


# Comments are removed.
sub recordBodyAsOneLine{
    my $this = shift;
    my $recordLinesRef = $this->get("_recordLines");
    return ""  unless( @$recordLinesRef );
    my $joinedLines = "";
    for( my $lineNo = 1; $lineNo < @$recordLinesRef; ++$lineNo ){
	my $line = @$recordLinesRef[$lineNo];
	next if $this->_isComment( $line );
	chomp $line;
	$joinedLines .= $line;
    }
    return $joinedLines;
}


sub recordCommentLinesRef{
    my $this = shift;
    my @recordBody = $this->recordBody();
    my $retVal = [];
    for (@recordBody){
	if( $this->_isComment($_) ){
	    push @$retVal, $_;
	}
    }
    return $retVal;
}


# Get a fasta format record and returns the lines it includes.
# Returns () if at end of stream.
#
# Use this method to iterate over fasta records.
sub getRecord{
    my $this = shift;
    my $istream = shift || *STDIN;

    return ()  if( $this->get("_atEndOfStream") );

#   Initialize at start of record.
    my $recordLinesRef = $this->get("_recordLines");
    @$recordLinesRef = ();
    my $recordWarningsRef = $this->get("_recordWarnings");
    @$recordWarningsRef = ();

    push @$recordLinesRef, $this->{_prevRecordLine} if( $this->get("_prevRecordLine") );
    my $headLineAlreadySeen = $this->get("_prevRecordLine");

    while( <$istream> ){
	my $curLine = $_;
	if( $this->_isComment($curLine) ){
	    unless( $headLineAlreadySeen ){
 		push @{$this->get("_streamHeaderLines")}, $curLine;
		next;
	    }
	}
	if( $curLine =~ /^>/ ){
	    $headLineAlreadySeen = 1;
	    ( $curLine =~ /^>\s*\S+/ ) || $this->_warn( 1, "Warning: A record with no record id encountered." );
	    if( $this->get("strip-white-space-before-id") ){
		$curLine =~ s/^>\s+/>/;
	    }
	    if(  !($this->get("_prevRecordLine") )  ){ # if first record in the stream
		push @$recordLinesRef, $curLine;
		$this->_set( "_prevRecordLine", $curLine );
		next;
	    }
	    # else this is the first line of the next record in the stream.
	    $this->_set( "_prevRecordLine", $curLine );

	    $this->_handleEndOfRecord( $recordLinesRef );
	    return @$recordLinesRef;
	}
	unless( $headLineAlreadySeen ){
	    chomp $curLine;
	    die( "Error: first line of what was expected to be a record did not start with \">\"\n"
		 . "input line: \"$curLine\"\n" );
	}
	$curLine =~ s|\s||g; # get rid of white space in sequences.
	push @$recordLinesRef, $curLine;
    }
    $this->_set( "_atEndOfStream", 1 );
    $this->_handleEndOfRecord( $recordLinesRef );
    return @$recordLinesRef;
}

sub printStreamHeader{
    my $this = shift;
    my $ostream = shift || *STDOUT;

    print( "Reader.pm: first line of _streamHeaderLines: ",
	   ${$this->{_streamHeaderLines}}, "\n" );
    
    print $ostream $_  for (@{$this->{_streamHeaderLines}});
}

# uses single hash as an argument. See %defaults below.
#
# Can print any sequence by passing it and optional record comment lines
#
# comment lines are printed after the head. Not reformatted.
# Default behavior is to print the current record.
sub printRecordFormated{
    my $this = shift;
    my %default = ( ostream => *STDOUT,
		     bodyLineLength => 70, 
		     recordBodyAsOneLine => $this->recordBodyAsOneLine(),
		     recordCommentLinesRef => $this->recordCommentLinesRef() );
    my %userOpt = @_;
    for my $userOpt (keys %userOpt){
	unless(  defined( $default{$userOpt} )  ){
	    Carp::croak( "Catastrophic Error, printRecordFormated: unknow option \"$userOpt\"" );
	  }
    }
    my %opts = ( %default, %userOpt );
    my $bodyLineLength = $opts{bodyLineLength};
    my $ostream = $opts{ostream};
    my $recordBodyAsOneLine = $opts{recordBodyAsOneLine};
    print $ostream $this->recordHead(), "\n";
    unless( $this->get("strip-comments") ){
	print $ostream "$_" for(@{$opts{recordCommentLinesRef}});
    }
    for( my $startPos = 0; $startPos < length( $recordBodyAsOneLine ); $startPos += $bodyLineLength ){
	print $ostream substr( $recordBodyAsOneLine, $startPos, $bodyLineLength );
	print $ostream "\n";
    }
}

sub printRecord{
    my $this = shift;
    my $ostream = shift || *STDOUT;
    my $recordLinesRef = $this->get("_recordLines");

    for my $line (@$recordLinesRef){
	unless( $this->get("strip-comments") && $this->_isComment($line) ){
	    print $ostream $line;
	}
    }
}


sub get{
    my ($this, $key) = @_;
    defined( $this->{$key} ) || Carp::croak( "Error: no data defined for key: \"$key\"" );
    return $this->{$key}
}


sub _set{
    my ($this, $key, $value) = @_;
    defined( $this->{$key} ) || Carp::croak( "Error: no data defined for key: \"$key\"" );
    $this->{$key} = $value;
}

sub _warn{
    my( $this, $level, $message ) = @_;
    if( $this->get("warn-level") >= $level ){
	push @{$this->get("_recordWarnings")}, "$message";
    }
}

sub _isComment{
    my ($this, $line) = @_;
    return 0  if( $this->get("no-comment") );
    return ( $line =~ /$this->{'comment-regex'}/ );
}


sub _handleEndOfRecord{
    my( $this, $recordLinesRef ) = @_;
    my $lastLineRef = \@$recordLinesRef[$#$recordLinesRef];

    my $id = $this->recordId();
    my @recordBody = $this->recordBody();
    unless( @recordBody ){
	$this->_warn( 0, "Warning: empty record, with id: $id, was encoutered." );
    }
    if( $this->get("add-asterisk") ){
	if( @recordBody ){
	    $$lastLineRef =~ s/([^*\n])$/*/;
	}else{ # empty body.
	    $this->_warn( 0, "Warning: empty record, with id: $id, was encoutered." );
	    push @$recordLinesRef, "*\n";
	}
	return;
    }

    # $this->get("add-asterisk") is false.
    @$recordLinesRef || return;  # if empty body do nothing.

    if( $this->get("strip-asterisk") ){
	$lastLineRef =~ s/\*$//;
    }
}


############################################################3
# Begin Class Methods
############################################################3

# class method. For calling from outside.
sub prettyPrintHashValues{
    shift; # get rid of the leading class name arg.
    return _prettyPrintHashValues( @_ );
}

# class method.
sub _prettyPrintHashValues{
    my $ostream = shift;
    my $hashRef = shift;
    my @keyNamesToPrint = @_;

    my $keyNameFieldLength = _maxStringLengthInList( @keyNamesToPrint );
    for my $keyName (@keyNamesToPrint){
	unless(  defined( $hashRef->{$keyName} )  ){
	    Carp::croak( "Catastrophic Error: no definition for keyName\"$keyName\"\n May merit submitting a bug report." );
	}
	my $valueNameForDisplay = $hashRef->{$keyName};
	$valueNameForDisplay =~ s/\*main:://;
	my $formatString = "  %-" . $keyNameFieldLength . "s %s\n";
	printf $ostream $formatString, $keyName, $valueNameForDisplay;
    }
}

# class method.
sub _maxStringLengthInList{
    my @l = @_;
    my $maxLength = 0;
    for my $s (@l){
	my $len = length $s;
	if( $len > $maxLength ){ $maxLength = $len; }
    }
    return $maxLength;
}


1;
