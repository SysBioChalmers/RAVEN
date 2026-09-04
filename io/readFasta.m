function fastaStruct=readFasta(fileName)
% readFasta  Read sequences from a FASTA file.
%
% Reads a FASTA file into a structure array, providing the same Header and
% Sequence fields as the Bioinformatics Toolbox fastaread, but using only
% base MATLAB so that no toolbox is required.
%
% Parameters
% ----------
% fileName : char
%     path to the FASTA file to read.
%
% Returns
% -------
% fastaStruct : struct
%     structure array with one element per sequence, with the fields:
%
%     - Header : the header line, without the leading ">"
%     - Sequence : the sequence with all line breaks removed
%
% Examples
% --------
%     fastaStruct = readFasta('proteins.fa');
%
% See also
% --------
% writeFasta

fileName=char(fileName);
fid=fopen(fileName,'r');
if fid==-1
    error('Cannot open FASTA file "%s".',fileName);
end

fastaStruct=struct('Header',{},'Sequence',{});
header='';
seq='';
haveRecord=false;
tline=fgetl(fid);
while ischar(tline)
    tline=strtrim(tline);
    if ~isempty(tline) && tline(1)=='>'
        if haveRecord
            fastaStruct(end+1).Header=header; %#ok<AGROW>
            fastaStruct(end).Sequence=seq;
        end
        header=strtrim(tline(2:end));
        seq='';
        haveRecord=true;
    elseif ~isempty(tline)
        seq=[seq tline]; %#ok<AGROW>
    end
    tline=fgetl(fid);
end
if haveRecord
    fastaStruct(end+1).Header=header;
    fastaStruct(end).Sequence=seq;
end
fclose(fid);
end
