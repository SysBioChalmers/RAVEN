function equationStrings=constructEquations(model,rxns,useComps,sortRevRxns,sortMetNames,useMetID,useFormula,useRevField)
% constructEquations
%   Construct equation strings for reactions
%
% Input:
%   model             a model structure
%   rxns              either a cell array of reaction IDs, a logical vector
%                     with the same number of elements as reactions in the
%                     model, or a vector of reaction indexes (optional, default
%                     model.rxns)
%   useComps          include the compartment of each metabolite (optional,
%                     default true)
%   sortRevRxns       sort reversible reactions so that the metabolite that
%                     is first in the lexiographic order is a reactant
%                     (optional, default false)
%   sortMetNames      sort the metabolite names in the equation. Uses
%                     compartment even if useComps is false (optional, default
%                     false)
%   useMetID          use metabolite ID in generated equations (optional,
%                     default false)
%   useFormula        use metabolite formula in generated equations (optional,
%                     default false)
%   useRevField       use the model.rev field to indicate reaction
%                     reversibility, alternatively this is determined from
%                     the model.ub and model.lb fields (optional, default true)
%
% Output:
%   equationStrings   a cell array with equations
%
% If useRevField is false, then reactions should be organized in their
% forward direction (e.g. ub = 1000 and lb = -1000/0) for the 
% reversibility to be correctly determined.
%
% Usage: equationStrings = constructEquations(model, rxns, useComps,...
%           sortRevRxns, sortMetNames, useMetID, useFormula, useRevField)

if nargin<2 || isempty(rxns)
    rxns=model.rxns;
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
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
if nargin<8
    useRevField=true;
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
    if useRevField == true
        isrev        = model.rev(Rindexes(i))==1;
    else
        isrev        = model.lb(Rindexes(i))<0 & model.ub(Rindexes(i))>0;
    end
    
    %Construct equation:
    equationStrings{i} = buildEquation(mets,stoichCoeffs,isrev);
end
end
