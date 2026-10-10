function writeFasta(fileName,fastaStruct,lineWidth)
% writeFasta  Write sequences to a FASTA file.
%
% Writes a structure array of sequences to a FASTA file, replacing the
% Bioinformatics Toolbox fastawrite with base MATLAB so that no toolbox is
% required. An existing file is overwritten.
%
% Parameters
% ----------
% fileName : char
%     path to the FASTA file to write.
% fastaStruct : struct
%     structure array with one element per sequence, with the fields:
%
%     - Header : the header line (written after a leading ">")
%     - Sequence : the sequence
% lineWidth : double
%     number of residues per sequence line (default 80).
%
% Examples
% --------
%     writeFasta('proteins.fa', fastaStruct);
%
% See also
% --------
% readFasta

if nargin<3 || isempty(lineWidth)
    lineWidth=80;
end
fileName=char(fileName);
fid=fopen(fileName,'w');
if fid==-1
    error('Cannot open FASTA file "%s" for writing.',fileName);
end

for i=1:numel(fastaStruct)
    fprintf(fid,'>%s\n',fastaStruct(i).Header);
    seq=fastaStruct(i).Sequence;
    for p=1:lineWidth:numel(seq)
        fprintf(fid,'%s\n',seq(p:min(p+lineWidth-1,numel(seq))));
    end
end
fclose(fid);
end
