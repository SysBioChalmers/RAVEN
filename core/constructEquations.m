function equationStrings=constructEquations(model,rxns,useComps,sortRevRxns,sortMetNames,useMetID)
% constructEquations
%   Construct equation strings for reactions
%
%   model             a model structure
%   rxns              either a cell array of reaction IDs, a logical vector with the
%                     same number of elements as reactions in the model, or a vector
%                     of reaction indexes (opt, default model.rxns)
%   useComps          include the compartment of each metabolite (opt, default true)
%   sortRevRxns       sort reversible reactions so that the metabolite that is first in
%                     the lexiographic order is a reactant (opt, default
%                     false)
%   sortMetNames      sort the metabolite names in the equation. Uses
%                     compartment even if useComps is false (opt, default
%                     false)
%   useMetID          use metabolite ID in generated equations, otherwise metNames are
%                     used (opt, default false)
%
%   equationStrings   a cell array with equations
%
%    Usage: equationStrings=constructEquations(model,rxns,useComps,...
%           sortRevRxns,sortMetNames,useMetID)
%
%   Hao Wang, 2017-05-15
%   Benjamin Sanchez, 2018-08-22
%

if nargin<2
    rxns=model.rxns;
end
if nargin<3
    useComps=true;
end
if nargin<4
    sortRevRxns=false;
end
if nargin<5
    sortMetNames=false;
end
if nargin<6
    useMetID=false;
end
if isempty(rxns) && nargin>2
    rxns=model.rxns;
end

%Sort reversible equations
if sortRevRxns==true
    model=sortModel(model);
end

%Sort metabolite names, including compartment
if sortMetNames==true
    model=sortModel(model,false,true);
end

indexes=getIndexes(model,rxns,'rxns');

equationStrings=cell(numel(indexes),1);

for i=1:numel(indexes)
    %Define metabolites by id or name, and with or without compartment:
    if useMetID==true
        mets = model.mets(model.S(:,indexes(i))~=0);
    else
        mets = model.metNames(model.S(:,indexes(i))~=0);
    end
    if useComps==true
        comps = model.comps{model.metComps(model.S(:,indexes(i))~=0)};
        mets  = strcat(mets,'[',comps,']');
    end
    %Define stoich coeffs and reversibility:
    stoichCoeffs = model.S(model.S(:,indexes(i))~=0,i);
    isrev        = model.rev(indexes(i))==1;
    
    %Construct equation:
    equationStrings{i} = buildEquation(mets,stoichCoeffs,isrev);
end

end
