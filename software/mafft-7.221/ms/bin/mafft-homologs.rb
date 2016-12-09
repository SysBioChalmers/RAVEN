#!/usr/bin/env ruby

localdb = "sp"        
# database name from which homologues are collected 
# by locally installed blast. Leave this if you do 
# not use the '-l' option.

mafftpath = "/usr/local/bin/mafft"   
# path of mafft. "/usr/local/bin/mafft"
# if mafft is in your command path, "mafft" is ok.

blastpath = "blastall"   
# path of blastall. 
# if blastall is in your command path, "blastall" is ok.

# mafft-homologs.rb  v. 2.1 aligns sequences together with homologues 
# automatically collected from SwissProt via NCBI BLAST.
#
# mafft > 5.58 is required
#
# Usage:
#   mafft-homologs.rb [options] input > output
# Options:
#   -a #      the number of collected sequences (default: 50)
#   -e #      threshold value (default: 1e-10)
#   -o "xxx"  options for mafft 
#             (default: " --op 1.53 --ep 0.123 --maxiterate 1000")
#   -l        locally carries out blast searches instead of NCBI blast
#             (requires locally installed blast and a database)
#   -f        outputs collected homologues also (default: off)
#   -w        entire sequences are subjected to BLAST search 
#             (default: well-aligned region only)

require 'getopts'
require 'tempfile'

# mktemp
GC.disable
temp_vf = Tempfile.new("_vf").path
temp_if = Tempfile.new("_if").path
temp_pf = Tempfile.new("_pf").path
temp_af = Tempfile.new("_af").path
temp_qf = Tempfile.new("_qf").path
temp_bf = Tempfile.new("_bf").path
temp_rid = Tempfile.new("_rid").path
temp_res = Tempfile.new("_res").path


system( mafftpath + " --help > #{temp_vf} 2>&1" )
pfp = File.open( "#{temp_vf}", 'r' )
while pfp.gets
	break if $_ =~ /MAFFT v/
end
pfp.close
if( $_ ) then
	mafftversion = sub( /^\D*/, "" ).split(" ").slice(0).strip.to_s
else
	mafftversion = "0"
end
if( mafftversion < "5.58" ) then
	puts ""
	puts "======================================================"
	puts "Install new mafft (v. >= 5.58)"
	puts "======================================================"
	puts ""
	exit
end

srand ( 0 )

def readfasta( fp, name, seq )
	nseq = 0
	tmpseq = ""
	while fp.gets
		if $_ =~ /^>/ then
			name.push( $_.sub(/>/,"").strip )
			seq.push( tmpseq ) if nseq > 0
			nseq += 1
			tmpseq = ""
		else
			tmpseq += $_.strip
		end
	end
	seq.push( tmpseq )
	return nseq
end

nadd = 50
eval = 1e-10
local = 0
fullout = 0
entiresearch = 0
corewin = 50
corethr = 0.3
mafftopt = " --op 1.53 --ep 0.123 --localpair --maxiterate 1000 --reorder "
if getopts( "s", "f", "w", "l", "h", "e:", "a:", "o:", "c:", "d:" ) == nil ||  ARGV.length == 0 || $OPT_h then
	puts "Usage: #{$0} [-h -l -e# -a# -o\"[options for mafft]\"] input_file"
	exit
end

if $OPT_c then
	corewin = $OPT_c.to_i
end
if $OPT_d then
	corethr = $OPT_d.to_f
end
if $OPT_w
	entiresearch = 1
end
if $OPT_f
	fullout = 1
end
if $OPT_s
	fullout = 0
end
if $OPT_l
	local = 1
end
if $OPT_e then
	eval = $OPT_e.to_f
end
if $OPT_a then
	nadd = $OPT_a.to_i
end
if $OPT_o then
	mafftopt += " " + $OPT_o + " "
end

system "cat " + ARGV.to_s + " > #{temp_if}"
ar = mafftopt.split(" ")
nar = ar.length
for i in 0..(nar-1)
	if ar[i] == "--seed" then
		system "cat #{ar[i+1]} >> #{temp_if}"
	end
end

nseq = 0
ifp = File.open( "#{temp_if}", 'r' )
	while ifp.gets
		nseq += 1 if $_ =~ /^>/
	end
ifp.close

if nseq >= 100 then
	STDERR.puts "The number of input sequences must be <100."
	exit
elsif nseq == 1 then
	system( "cp #{temp_if}"  + " #{temp_pf}" )
else
	STDERR.puts "Performing preliminary alignment .. "
	if entiresearch == 1 then
#		system( mafftpath + " --maxiterate 1000 --localpair #{temp_if} > #{temp_pf}" )
		system( mafftpath + " --maxiterate 0 --retree 2 #{temp_if} > #{temp_pf}" )
	else
		system( mafftpath + " --maxiterate 1000 --localpair --core --coreext --corethr #{corethr.to_s} --corewin #{corewin.to_s} #{temp_if} > #{temp_pf}" )
	end
end

pfp = File.open( "#{temp_pf}", 'r' )
inname = []
inseq = []
slen = []
act = []
nin = 0
nin = readfasta( pfp, inname, inseq )
for i in 0..(nin-1)
	slen.push( inseq[i].gsub(/-/,"").length )
	act.push( 1 )
end
pfp.close

pfp = File.open( "#{temp_if}", 'r' )
orname = []
orseq = []
nin = 0
nin = readfasta( pfp, orname, orseq )
pfp.close

allen = inseq[0].length
for i in 0..(nin-2)
	for j in (i+1)..(nin-1)
		next if act[i] == 0
		next if act[j] == 0
		pid = 0.0
		total = 0
		for a in 0..(allen-1)
			next if inseq[i][a,1] == "-" || inseq[j][a,1] == "-"
			total += 1
			pid += 1.0 if inseq[i][a,1] == inseq[j][a,1]
		end
		pid /= total
#		puts "#{i.to_s}, #{j.to_s}, #{pid.to_s}"
		if pid > 0.5 then
			if slen[i] < slen[j]
				act[i] = 0 
			else
				act[j] = 0 
			end
		end
	end
end
#p act


afp = File.open( "#{temp_af}", 'w' )

STDERR.puts "Searching .. \n"
ids = []
add = []
sco = []
for i in 0..(nin-1)
	inseq[i].gsub!(/-/,"")
	afp.puts ">" + orname[i]
	afp.puts orseq[i]

#	afp.puts ">" + inname[i]
#	afp.puts inseq[i]

	STDERR.puts "Query (#{i+1}/#{nin})\n" + inname[i]
	if act[i] == 0 then
		STDERR.puts "Skip.\n\n"
		next 
	end

	if local == 0 then
		command = "lynx -source 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?QUERY=" + inseq[i] + "&DATABASE=swissprot&HITLIST_SIZE=" + nadd.to_s + "&FILTER=L&EXPECT='" + eval.to_s + "'&FORMAT_TYPE=TEXT&PROGRAM=blastp&SERVICE=plain&NCBI_GI=on&PAGE=Proteins&CMD=Put' > #{temp_rid}"
		system command
	
		ridp = File.open( "#{temp_rid}", 'r' )
		while ridp.gets
			break if $_ =~ / RID = (.*)/
		end
		ridp.close
		rid = $1.strip
		STDERR.puts "Submitted to NCBI. rid = " + rid
	
		STDERR.printf "Waiting "
		while 1 
			STDERR.printf "."
			sleep 10
			command = "lynx -source 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RID=" + rid + "&DESCRIPTIONS=500&ALIGNMENTS=" + nadd.to_s + "&ALIGNMENT_TYPE=Pairwise&OVERVIEW=no&CMD=Get&FORMAT_TYPE=XML' > #{temp_res}"
			system command
			resp = File.open( "#{temp_res}", 'r' )
#			resp.gets
#			if $_ =~ /WAITING/ then
#				resp.close
#				next
#			end
			while( resp.gets )
				break if $_ =~ /QBlastInfoBegin/
			end
			resp.gets
			if $_ =~ /WAITING/ then
				resp.close
				next
			else
				resp.close
				break
			end
		end
	else
#		puts "Not supported"
#		exit
		qfp = File.open( "#{temp_qf}", 'w' )
			qfp.puts "> "
			qfp.puts inseq[i]
		qfp.close
		command = blastpath + "  -p blastp  -e #{eval} -b 1000 -m 7 -i #{temp_qf} -d #{localdb} > #{temp_res}"
		system command
		resp = File.open( "#{temp_res}", 'r' )
	end
	STDERR.puts " Done.\n\n"

	resp = File.open( "#{temp_res}", 'r' )
	while 1
		while resp.gets
			break if $_ =~ /<Hit_id>(.*)<\/Hit_id>/ || $_ =~ /(<Iteration_stat>)/
		end
		id = $1
		break if $_ =~ /<Iteration_stat>/
#		p id
		while resp.gets
			break if $_ =~ /<Hsp_bit-score>(.*)<\/Hsp_bit-score>/
		end
		score = $1.to_f
#		p score

		known = ids.index( id )
		if known != nil then
			if sco[known] >= score then
				next
			else
				ids.delete_at( known )
				add.delete_at( known )
				sco.delete_at( known )
			end
		end
		while resp.gets
			break if $_ =~ /<Hsp_hseq>(.*)<\/Hsp_hseq>/
		end
#		break if $1 == nil
		target = $1.sub( /-/, "" ).sub( /U/, "X" )
#		p target
#		STDERR.puts "adding 1 seq"
		ids.push( id )
		sco.push( score )
		add.push( target )
	end
	resp.close
end

n = ids.length

outnum = 0
while n > 0 && outnum < nadd
	m = rand( n )
	afp.puts ">_addedbymaffte_" + ids[m]
	afp.puts add[m]
	ids.delete_at( m )
	add.delete_at( m )
	n -= 1
	outnum += 1
end
afp.close

STDERR.puts "Performing alignment .. "
system( mafftpath + mafftopt + " #{temp_af} > #{temp_bf}" )
STDERR.puts "done."

bfp = File.open( "#{temp_bf}", 'r' )
outseq = []
outnam = []
readfasta( bfp, outnam, outseq )
bfp.close

outseq2 = []
outnam2 = []

len = outseq.length
for i in 0..(len-1)
#	p outnam[i]
	if fullout == 0 && outnam[i] =~ /_addedbymaffte_/ then
		next
	end
	outseq2.push( outseq[i] )
	outnam2.push( outnam[i].sub( /_addedbymaffte_/, "_ho_" ) )
end

nout = outseq2.length
len = outseq[0].length
p = len
while p>0
	p -= 1
    allgap = 1
    for j in 0..(nout-1)
		if outseq2[j][p,1] != "-" then
			allgap = 0
			break
		end
    end
    if allgap == 1 then
        for j in 0..(nout-1)
            outseq2[j][p,1] = ""
        end
    end
end
for i in 0..(nout-1)
	puts ">" + outnam2[i]
	puts outseq2[i].gsub( /.{1,60}/, "\\0\n" )
end


system( "rm -rf #{temp_if} #{temp_vf} #{temp_af} #{temp_bf} #{temp_pf} #{temp_qf} #{temp_res} #{temp_rid}" )
