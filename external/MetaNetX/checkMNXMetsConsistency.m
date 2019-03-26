function result=checkMNXMetsConsistency(metMetaNetXID,MNXmets)
% checkMNXMetsConsistency
%   checkMNXMetsConsistency checks consistency in chemical formula and
%   charge when multiple MNX metids are associated to the same metabolite.
%
%   metMetaNetXID    cell array derived from model structure, with MNX
%               metabolite ids.
%   MNXmets     MNXref structure, reconstructed by buildMNXref. (opt, if
%               not specified, buildMNXref is called)
%
%   result      vector of scores with length, detailing:
%               0   no MNXid annotated for this metabolite
%               1   one MNXid annotated for this metabolite
%               2   multiple MNXids anotated for this metabolite, all with
%                   the same chemical formula and charge
%               -1  multiple MNXids anotated for this metabolite, with
%                   different chemical formula and/or charge
%
%   Usage: result=checkMNXMetsConsistency(metMetaNetXID,MNXmets)
%
%   Eduard Kerkhoven, 2018-07-16

if nargin<2
    MNXmets=buildMNXref('mets');
end

% Sort MNXMets structure to speedup later processes
[MNXmets.metMetaNetXID, I]=sort(MNXmets.metMetaNetXID);
fields=fieldnames(MNXmets);
for i=1:length(fields)
    if isequal(size(eval(strcat('MNXmets.',fields{i})),1),length(I))
        eval(strcat('MNXmets.',fields{i},'=MNXmets.',fields{i},'(I);'));
    end
end

% Initilize output cell array
result=zeros(size(metMetaNetXID,1),1);

% no MNX metID associated: score 0; one MNX metID associated: score 1
numIds=cellfun(@numel,metMetaNetXID);
%numIds(numIds>1) = 1;
%numIds = sum(numIds,2);
% result(find(numIds==0))=0; % already 0
result(numIds==1)=1;

% loop through multiple MNX associations
numIds=find(numIds>1);
for i=1:length(numIds)
    query=metMetaNetXID{numIds(i),:};
	[hit, j]=ismember(query,MNXmets.mets);
    if all(hit)
        if numel(hit)>1
            charges=num2cell(MNXmets.metCharges(j));
            if isequal(MNXmets.metFormulas{j}) && isequal(charges{:})
                result(numIds(i))=2;
            else
                result(numIds(i))=-1;
            end
		end
	else
		error('The model contains unknown MNX identifiers.');
	end
end
end