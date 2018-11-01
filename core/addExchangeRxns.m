function [model, addedRxns]=addExchangeRxns(model,reactionType,mets)
% addExchangeRxns
%   Adds exchange reactions for some metabolites
%
%   model           a model structure
%   reactionType    the type of reactions to add
%                   'in'    input reactions
%                   'out'   output reactions
%                   'both'  reversible input/output reactions. Positive
%                   direction corresponds to output
%   mets            either a cell array of metabolite IDs, a logical vector
%                   with the same number of elements as metabolites in the model,
%                   or a vector of indexes to add for (opt, default model.mets)
%
%   model           updated model structure
%   addedRxns       ids of the added reactions
%
%   This is a faster version than addRxns when adding exchange reactions.
%   New reactions are named "metName exchange (OUT/IN/BOTH)" while reaction
%   ids are formatted as "EXC_OUT/IN/BOTH_METID".
%
%   Usage: [model, addedRxns]=addExchangeRxns(model,reactionType,mets)
%
%   Simonas Marcisauskas, 2018-03-17
%

if nargin<3
    mets=model.mets;
end
reactionType=upper(reactionType);
J=getIndexes(model,mets,'mets',false);
mets=model.mets(J);

%Production is positive for OUT and BOTH
if strcmp(reactionType,'IN')
    I=speye(numel(model.mets));
else
    I=speye(numel(model.mets))*-1;
end
I=I(:,J);

%Add an I matrix which corresponds to production of all metabolites
model.S=[model.S I];
if strcmp(reactionType,'BOTH')
    model.lb=[model.lb;ones(numel(J),1)*-1000];
    model.rev=[model.rev;ones(numel(J),1)];
else
    model.lb=[model.lb;zeros(numel(J),1)];
    model.rev=[model.rev;zeros(numel(J),1)];
end
model.ub=[model.ub;ones(numel(J),1)*1000];
model.c=[model.c;zeros(numel(J),1)];

filler=cell(numel(J),1);
filler(:)={''};
addedRxns=strcat('EXC_',reactionType,'_',mets);
addedRxnNames=strcat(model.metNames(J),' exchange (',reactionType,')');
model.rxns=[model.rxns;addedRxns];
model.rxnNames=[model.rxnNames;addedRxnNames];

if isfield(model,'eccodes')
    model.eccodes=[model.eccodes;filler];
end
if isfield(model,'subSystems')
    model.subSystems=[model.subSystems;filler];
end
if isfield(model,'grRules')
    model.grRules=[model.grRules;filler];
end
if isfield(model,'rxnFrom')
    model.rxnFrom=[model.rxnFrom;filler];
end
if isfield(model,'rxnMiriams')
    model.rxnMiriams=[model.rxnMiriams;filler];
end
if isfield(model,'rxnGeneMat')
    model.rxnGeneMat=[model.rxnGeneMat;sparse(numel(J),numel(model.genes))];
end
if isfield(model,'rxnComps')
    model.rxnComps=[model.rxnComps;ones(numel(J),1)];
    fprintf('NOTE: The exchange reactions are assigned to the first compartment\n');
end
if isfield(model,'rxnNotes')
    model.rxnNotes=[model.rxnNotes;filler];
end
if isfield(model,'rxnReferences')
    model.rxnReferences=[model.rxnReferences;filler];
end
if isfield(model,'rxnConfidenceScores')
    model.rxnConfidenceScores=[model.rxnConfidenceScores;NaN(numel(J),1)];
end
end
