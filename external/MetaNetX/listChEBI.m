function listChEBI()
% listSecondaryChEBI
%   Makes a list of primary ChEBI identifiers, to filter out secondary and
%   redundant ChEBI IDs.
%
% Usage: listChEBI(path)
%
% Eduard Kerkhoven, 2018-09-24

[ST, I]=dbstack('-completenames');
filePath=fileparts(fileparts(fileparts(ST(I).file)));
filePath=fullfile(filePath,'external','metanetx');

if ~exist(fullfile(filePath,'chebi_lite.obo'), 'file')
    error(['File chebi_lite.obo cannot be found please download from\n'...
        'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo\n'...
        'and store in %s'],filePath);
end

fid=fopen('chebi_lite.obo');
chebi=fscanf(fid,'%s');
fclose(fid);
chebi2=chebi;
chebi=strsplit(chebi,'[Term]');
chebi(1)=[];
chebi=regexprep(chebi,'^id:CHEBI:','');
chebi=regexprep(chebi,'^(\d+).*','$1');
chebi=sort(str2double(chebi));

fid=fopen(fullfile(filePath,'chebi.dat'),'w');
fprintf(fid,'%i\n',chebi);
fclose(fid);
end
