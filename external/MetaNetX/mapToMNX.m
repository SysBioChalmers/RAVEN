function [MNXfields, MetMNXmatch]=mapToMNX(model,rxns,MNXref,keepOneMetMNX)
% mapToMNX
%   Maps a model to MetaNetX as much as possible, 
%
%   model           model structure
%   rxns            logical whether reactions should be mapped to MetaNetX.
%                   'false' matches only metabolites, 'true' matches both
%                   metabolites and reactions (opt, default 'true')
%   MNXref          MNXref structure, reconstructed by buildMNXmodel. (opt,
%                   if not specified, buildMNXmodel is called)
%   keepOneMetMNX   logical, whether metabolites should only retain one
%                   MNXid. If 'true', then metabolites with multiple MNXids
%                   will be checked for identical chemical formula and
%                   charge. If identical, the earliest MNXid is kept, while
%                   alternative situations are detailed in MetMNXmatch.
%                   (opt, default 'true' if MetMNXmatch is specified as
%                   output)
%   
%   MNXfields       structure containing metMNXID and rxnMNXID, depending
%                   'rxns' parameter setting
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
        MNXref=buildMNXmodel('both');
    else
        MNXref=buildMNXmodel('mets');
    end
end
if ~exist('keepOneMetMNX','var') && nargout<2
    keepOneMetMNX=false;
else
    keepOneMetMNX=true;
end

if ~isfield(model,'rules') % If RAVEN model, first convert to Cobra    
    fprintf('Extract annotations by converting to COBRA model.\n');
    model=ravenCobraWrapper(model);
end

%Initial met to MNX mapping
model=mapModelMets(model,MNXref);

%Map rxns to MNX via metabolites, and filter metMNXIDs with these
%rxnMNXIDs.
if rxns
    rxnMaps=mapRxnsViaMets(model,MNXref);
    model.rxnMNXID=rxnMaps.rxnMNXID;
    [model,~] = filterMetMNXIDsViaRxns(model,MNXref,true,true);
end

if keepOneMetMNX
    result=checkMNXMetsConsistency(model.metMNXID,MNXref);
    newMNXID=strings(size(model.mets));
    oneMatch=find(result==1);
    newMNXID(oneMatch)=model.metMNXID(oneMatch,1);
    twoMatch=find(result==2);
    for i=1:length(twoMatch)
        MNXid=model.metMNXID{twoMatch(i),:};
        MNXid=MNXid(~cellfun('isempty',MNXid));
        MNXid=sort(MNXid);
        newMNXID(twoMatch(i))=MNXid(1);
    end
    model.metMNXID=newMNXID;
    MetMNXmatch=result;
end
MNXfields.metMNXID=model.metMNXID;
if rxns
    MNXfields.rxnMNXID=model.rxnMNXID;
end
end


