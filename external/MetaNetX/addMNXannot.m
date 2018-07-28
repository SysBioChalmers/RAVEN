function newModel = addMNXannot(model,MNXfields,MNXref,fields)
% addMNXannot
%   checkMNXMetsConsistency checks consistency in chemical formula and
%   charge when multiple MNX metids are associated to the same metabolite.
%
%   model       a model structure, as used to run mapToMNX
%   MNXfields   output from mapToMNX
%   MNXref      MNXref structure, reconstructed by buildMNXref. (opt,
%               if not specified, buildMNXref is called)
%   fields      cell array specifying what annotation fields from MNXref
%               should be appended to model (opt, default appends most
%               relevant fields, as detailed below)
%
%   newModel    model with new annotations
%
%   Note:   Default MNXfields all include all annotation fields, except for:
%           (i)     metSABIORKID and metREACTOMEID, as these are internal
%                   identifiers of their respective databases, and not
%                   identifiers.org namespaces
%           (ii)    metEnviPathID due to limited usefulness for most GEMs
%           (iii)   metHMDBID as it is only relevant for human models
%           (iv)    rxnREACTOMEID as all IDs are species specific.
%           Default MNXfields does include rxnSEEDID, even while it is not
%           an identifiers.org namespace.
%
%           Note that rxnRheaID contains all four identifiers, representing
%           all possible reaction reversibilities.   
%
%   Usage: result=addMNXannot(model,MNXfields,fields)
%
%   Eduard Kerkhoven, 2018-07-27

if nargin<4
    fields={'rxnBIGGID','rxnKEGGID','rxnMetaCycID','rxnRheaID',...
        'rxnSABIORKID','rxnSEEDID','rxnMetaNetXID','metBIGGID','metChEBIID',...
        'metKEGGID','metLIPIDMAPSID','metMetaCycID','metSEEDID',...
        'metSLMID','metMetaNetXID'};
end

if ischar(fields)
    fields={fields};
end

if ~exist('MNXref','var')
    MNXref=buildMNXref('both');
end

%check if metMetaNetXID or rxnMetaNetXID should be added, as they can be appended
idx=find(~cellfun(@isempty, regexp(fields,'(met|rxn)MNXID')));
for i=1:length(idx)
    j=fields{idx(i)};
    model.(j)=MNXfields.(j);
end

% reactions
fieldIdx=find(~cellfun(@isempty, regexp(fields,'^rxn.*'))); %which fields are rxn

MNXID=flattenCell(flattenCell(MNXfields.rxnMetaNetXID,true));
[Lia, Locb]=ismember(MNXID,MNXref.rxns);

for j=1:length(fieldIdx)
    annotStruct=strings(numel(model.rxns),1);
    for i=1:size(MNXID,1)
        idx=Locb(i,:);
        if ~sum(idx)==0
            idx=idx(Lia(i,:));
            DBID=MNXref.(fields{fieldIdx(j)})(idx);
            DBID=DBID(~cellfun(@isempty,DBID));
            if ~isempty(DBID)
                DBID=flattenCell(DBID);
                DBID=DBID(~cellfun(@isempty,DBID));
                annotStruct{i}=strjoin(DBID,';');
            end
        end
    end
    model.(fields{fieldIdx(j)})=annotStruct;
end


fieldIdx=find(~cellfun(@isempty, regexp(fields,'^met.*'))); %which fields are met

MNXID=flattenCell(flattenCell(MNXfields.metMetaNetXID,true));
[Lia, Locb]=ismember(MNXID,MNXref.mets);

for j=1:length(fieldIdx)
    annotStruct=strings(numel(model.rxns),1);
    for i=1:size(MNXID,1)
        idx=Locb(i,:);
        if ~sum(idx)==0
            idx=idx(Lia(i,:));
            DBID=MNXref.(fields{fieldIdx(j)})(idx);
            DBID=DBID(~cellfun(@isempty,DBID));
            if ~isempty(DBID)
                DBID=flattenCell(DBID);
                DBID=DBID(~cellfun(@isempty,DBID));
                annotStruct{i}=strjoin(DBID,';');
            end
        end
    end
    model.(fields{fieldIdx(j)})=annotStruct;
end
newModel=model;
end
