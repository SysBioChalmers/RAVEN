function transportRxns=getTransportRxns(model)
% getTransportRxns
%   Retrieves the transport reactions from a model
%
%   model           a model structure
%
%   transportRxns   logical array with true if the corresponding 
%                   reaction is a transport reaction
%
%   Transport reactions are defined as reactions involving (at least) one
%   metabolite name in more than one compartment.
%
% Usage: transportRxns=getTransportRxns(model)

transportRxns=false(numel(model.rxns),1);

for i=1:numel(model.rxns)
    %Get the involved metabolites in each reaction
    mets=model.metNames(model.S(:,i)~=0);
    transportRxns(i)=numel(mets)~=numel(unique(mets));
end
end
