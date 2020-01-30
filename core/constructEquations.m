function equationStrings=constructEquations(model,rxns,useComps,sortRevRxns,sortMetNames,useMetID,useFormula)
% constructEquations
%   Construct equation strings for reactions
%
%   Input:
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
%   useMetID          use metabolite ID in generated equations (opt,
%                     default false)
%   useFormula        use metabolite formula in generated equations (opt,
%                     default false)
%
%   Outut:
%   equationStrings   a cell array with equations
%
%   NOTE: Reactions in a model should be organized in their forward direction
%   (e.g. ub = 1000 and lb = -1000/0) so that their equations can be correctly
%   constructed by this function.
%
%   Usage: equationStrings=constructEquations(model,rxns,useComps,...
%           sortRevRxns,sortMetNames,useMetID,useFormula)
%
%   Hao Wang, 2019-03-05
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
if nargin<7
    useFormula=false;
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

Rindexes=getIndexes(model,rxns,'rxns');

equationStrings=cell(numel(Rindexes),1);

for i=1:numel(Rindexes)
    Mindexes=find(model.S(:,Rindexes(i)));
    %Define metabolites by id, formula or name, and with or without compartment: 
    if useMetID==true && useFormula==false
        mets = model.mets(Mindexes);
    elseif useMetID==false && useFormula==true
        mets = strcat('[',model.metFormulas(Mindexes),']');
    elseif useMetID==true && useFormula==true
        error('Arguments useMetID and useFormula cannot be both TRUE!');
    else
        mets = model.metNames(Mindexes);
    end
    if useComps==true
        comps = model.comps(model.metComps(Mindexes));
        mets  = strcat(mets,'[',comps,']');
    end
    %Define stoich coeffs and reversibility:
    stoichCoeffs = model.S(Mindexes,Rindexes(i));
    isrev        = model.lb(Rindexes(i))<0 & model.ub(Rindexes(i))>0;
    
    %Construct equation:
    equationStrings{i} = buildEquation(mets,stoichCoeffs,isrev);
end

end
