function listSecondaryChEBI()
% listSecondaryChEBI
%   Makes a list of secondary ChEBI identifiers that can be left out of
%   MNXref database.
%
% Usage: listSecondaryChEBI(path)
%
% Eduard Kerkhoven, 2018-07-25
% 

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

chebi=strsplit(chebi,'alt_id:');
chebi(1)=[];
chebi=regexprep(chebi,'(CHEBI:\d+).*','$1');
chebi=regexprep(chebi,'CHEBI:','');
chebi=sort(str2double(chebi));

fid=fopen(fullfile(filePath,'chebiSecondary.dat'),'w');
fprintf(fid,'%i\n',chebi);
fclose(fid);
end
