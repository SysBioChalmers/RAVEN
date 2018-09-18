function [MNXfields, MetMNXmatch]=mapToMNX(model,rxns,MNXref,keepOneMetMNX)
% mapToMNX
%   Maps a model to MetaNetX as much as possible,
%
%   model           model structure
%   rxns            logical whether reactions should be mapped to MetaNetX.
%                   'false' matches only metabolites, 'true' matches both
%                   metabolites and reactions (opt, default 'true')
%   MNXref          MNXref structure, reconstructed by buildMNXref. (opt,
%                   if not specified, buildMNXref is called)
%   keepOneMetMNX   logical, whether metabolites should only retain one
%                   MNXid. If 'true', then metabolites with multiple MNXids
%                   will be checked for identical chemical formula and
%                   charge. If identical, the earliest MNXid is kept, while
%                   alternative situations are detailed in MetMNXmatch.
%                   (opt, default 'true' if MetMNXmatch is specified as
%                   output)
%
%   MNXfields       structure containing metMetaNetXID and/or rxnMetaNetXID,
%                   depending on the 'rxns' parameter setting
%   MetMNXmatch     vector
%
%   Usage: MNXfields=mapToMNX(model,rxns,MNXref)
%
%   Eduard Kerkhoven, 2018-07-16

if ~exist('rxns','var')
    rxns=true;
end

if ~exist('MNXref','var')
    if rxns
        MNXref=buildMNXref('both');
    else
        MNXref=buildMNXref('mets');
    end
end
if ~exist('keepOneMetMNX','var')
    if nargout>1
        keepOneMetMNX=false;
    else
        keepOneMetMNX=true;
    end
end

if ~any(isfield(model,{'metMiriams','rxnMiriams'})) % Convert Miriam annotations
    model=convertMiriams(model);
    fprintf('Converting Miriam annotations to COBRA-style...\n');
end

%Initial met to MNX mapping
model=mapModelMets(model,MNXref);

%Map rxns to MNX via metabolites, and filter metMetaNetXIDs with these
%rxnMetaNetXIDs.
if rxns
    rxnMaps=mapRxnsViaMets(model,MNXref);
    model.rxnMetaNetXID=rxnMaps.rxnMetaNetXID;
    [model,~] = filterMetMNXIDsViaRxns(model,MNXref,true,true);
end

if keepOneMetMNX
    result=checkMNXMetsConsistency(model.metMetaNetXID,MNXref);
    newMNXID=strings(size(model.mets));
    oneMatch=find(result==1);
    newMNXID(oneMatch)=model.metMetaNetXID(oneMatch,1);
    twoMatch=find(result==2);
    for i=1:length(twoMatch)
        MNXid=model.metMetaNetXID{twoMatch(i),:};
        MNXid=MNXid(~cellfun('isempty',MNXid));
        MNXid=sort(MNXid);
        newMNXID(twoMatch(i))=MNXid(1);
    end
    model.metMetaNetXID=newMNXID;
    MetMNXmatch=result;
end
MNXfields.metMetaNetXID=model.metMetaNetXID;
if rxns
    MNXfields.rxnMetaNetXID=model.rxnMetaNetXID;
end
end
