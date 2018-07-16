function result=checkMNXMetsConsistency(metMNXID,MNXmets)
% checkMNXMetsConsistency
%   checkMNXMetsConsistency checks consistency in chemical formula and
%   charge when multiple MNX metids are associated to the same metabolite.
%
%   metMNXID    cell array derived from model structure, with MNX
%               metabolite ids.
%   MNXmets     MNXref structure, reconstructed by buildMNXmodel. (opt, if
%               not specified, buildMNXmodel is called)
%
%   result      vector of scores with length, detailing:
%               0   no MNXid annotated for this metabolite
%               1   one MNXid annotated for this metabolite
%               2   multiple MNXids anotated for this metabolite, all with
%                   the same chemical formula and charge
%               -1  multiple MNXids anotated for this metabolite, with
%                   different chemical formula and/or charge
%
%   Usage: result=checkMNXMetsConsistency(metMNXID,MNXmets)
%
%   Eduard Kerkhoven, 2018-07-16

if nargin<2
    MNXmets=buildMNXmodel('mets');
end

% Sort MNXMets structure to speedup later processes
[MNXmets.metMNXID, I]=sort(MNXmets.metMNXID);
fields=fieldnames(MNXmets);
for i=1:length(fields)
    if isequal(size(eval(strcat('MNXmets.',fields{i})),1),length(I))
        eval(strcat('MNXmets.',fields{i},'=MNXmets.',fields{i},'(I);'));
    end
end

% Initilize output cell array
result=zeros(size(metMNXID,1),1);

% no MNX metID associated: score 0; one MNX metID associated: score 1
numIds=cellfun(@numel,metMNXID);
% result(find(numIds==0))=0; % already 0
result(find(numIds==1))=1;

% loop through multiple MNX associations
numIds=find(numIds>1);
for i=1:length(numIds)
    query=metMNXID{numIds(i),:};
	[hit, j]=ismember(query,MNXmets.mets);
    if all(hit)
        charges=num2cell(MNXmets.metCharges(j));
		if isequal(MNXmets.metFormulas{j}) && isequal(charges{:})
            result(numIds(i))=2;
        else
			result(numIds(i))=-1;
		end
	else
		error('The model contains unknown MNX identifiers.');
	end
end
end