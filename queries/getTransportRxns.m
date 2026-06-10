function transportRxns=getTransportRxns(model)
% getTransportRxns  Retrieve the transport reactions from a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Returns
% -------
% transportRxns : logical
%     logical array with true if the corresponding reaction is a transport
%     reaction.
%
% Examples
% --------
%     transportRxns = getTransportRxns(model);
%
% Notes
% -----
% Transport reactions are defined as reactions involving (at least) one
% metabolite name in more than one compartment.

transportRxns=false(numel(model.rxns),1);

for i=1:numel(model.rxns)
    %Get the involved metabolites in each reaction
    mets=model.metNames(model.S(:,i)~=0);
    transportRxns(i)=numel(mets)~=numel(unique(mets));
end
end
