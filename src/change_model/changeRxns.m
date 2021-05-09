function model=changeRxns(model,rxns,equations,eqnType,compartment,allowNewMets)
% changeRxns
%   Modifies the equations of reactions
%
%   model            a model structure
%   rxns             cell array with reaction ids
%   equations        cell array with equations. Alternatively, it can be a
%                    structure with the fields "mets" and "stoichCoeffs",
%                    in the same fashion as addRxns. E.g.:
%                    equations.mets = {{'met1','met2'},{'met1','met3'}}
%                    equations.stoichCoeffs = {[-1,+2],[-1,+1]}
%   eqnType          double describing how the equation string should be
%                    interpreted
%                    1 - The metabolites are matched to model.mets. New
%                        metabolites (if allowed) are added to
%                        "compartment"
%                    2 - The metabolites are matched to model.metNames and
%                        all metabolites are assigned to "compartment". Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are
%                        already used in the model then the numbering will
%                        start from the highest used integer+1
%                    3 - The metabolites are written as
%                        "metNames[compNames]". Only compartments in
%                        model.compNames are allowed. Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are
%                        already used in the model then the numbering will
%                        start from the highest used integer+1
%   compartment      a string with the compartment the metabolites should
%                    be placed in when using eqnType=2. Must match
%                    model.compNames (opt when eqnType=1 or eqnType=3)
%   allowNewMets     true if the function is allowed to add new
%                    metabolites. It is highly recommended to first add
%                    any new metabolites with addMets rather than
%                    automatically through this function. addMets supports
%                    more annotation of metabolites, allows for the use of
%                    exchange metabolites, and using it reduces the risk
%                    of parsing errors (opt, default false)
%
%   model            an updated model structure
%
%   NOTE: This function should be used with some care, since it doesn't
%   care about bounds on the reactions. Changing a irreversible reaction to
%   a reversible one (or the other way around) will only change the
%   model.rev field and not the model.lb/model.ub fields. The reaction will
%   therefore still be having the same reversibility because of the
%   bounds. Use setParams to change the bounds.
%
%   NOTE: When adding metabolites to a compartment where it previously
%   doesn't exist, the function will copy any available information from
%   the metabolite in another compartment.
%
%   Usage: model=changeRxns(model,rxns,equations,eqnType,compartment,allowNewMets)

if nargin<4 && isfield(equations,'stoichCoeffs')
    eqnType=1;
end

if nargin<5
    compartment=[];
end
if nargin<6
    allowNewMets=false;
end

if ischar(rxns)
    rxns={rxns};
end
if ischar(equations)
    equations={equations};
end

%Find the indexes of the reactions and throw an error if they aren't all
%found
[I, J]=ismember(rxns,model.rxns);
if ~all(I)
    EM='All reaction ids must exist in the model';
    dispEM(EM);
end

%The reactions are changed in the following manner. First create a
%rxns-structure by copying info from the model. Then remove the old
%reactions. Then add the updated ones using addRxns. Lastly, the model is
%reordered to match the original order. This is done like this to make use
%of the advanced parsing of equations that addRxns use.
rxnsToChange.rxns=rxns;
if isfield(equations,'mets') && isfield(equations,'stoichCoeffs')
    rxnsToChange.mets=equations.mets;
    rxnsToChange.stoichCoeffs=equations.stoichCoeffs;
else
    rxnsToChange.equations=equations;
end
if isfield(model,'rxnNames')
    rxnsToChange.rxnNames=model.rxnNames(J);
end
if isfield(model,'lb')
    rxnsToChange.lb=model.lb(J);
end
if isfield(model,'ub')
    rxnsToChange.ub=model.ub(J);
end
if isfield(model,'c')
    rxnsToChange.c=model.c(J);
end
if isfield(model,'eccodes')
    rxnsToChange.eccodes=model.eccodes(J);
end
if isfield(model,'subSystems')
    rxnsToChange.subSystems=model.subSystems(J);
end
if isfield(model,'rxnComps')
    rxnsToChange.rxnComps=model.rxnComps(J);
end
if isfield(model,'grRules')
    rxnsToChange.grRules=model.grRules(J);
end
if isfield(model,'rxnFrom')
    rxnsToChange.rxnFrom=model.rxnFrom(J);
end
if isfield(model,'rxnScores')
    rxnsToChange.rxnScores=model.rxnScores(J);
end
if isfield(model,'rxnMiriams')
    rxnsToChange.rxnMiriams=model.rxnMiriams(J);
end
if isfield(model,'rxnNotes')
    rxnsToChange.rxnNotes=model.rxnNotes(J);
end
if isfield(model,'rxnReferences')
    rxnsToChange.rxnReferences=model.rxnReferences(J);
end
if isfield(model,'rxnConfidenceScores')
    rxnsToChange.rxnConfidenceScores=model.rxnConfidenceScores(J);
end
if isfield(model,'pwys')
    rxnsToChange.pwys=model.pwys(J);
end

%Calculate the new order of reactions
order=ones(numel(model.rxns),1);
order(J)=0;
order=cumsum(order);
order(J)=order(end)+1:order(end)+numel(rxns);

%Remove the original reactions
model=removeReactions(model,rxns);

model=addRxns(model,rxnsToChange,eqnType,compartment,allowNewMets);
model=permuteModel(model,order,'rxns');
end
