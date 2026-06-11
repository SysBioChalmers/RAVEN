function equationStrings=constructEquations(model,varargin)
% constructEquations  Construct equation strings for reactions.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% rxns : cell or logical or double, optional
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of reaction
%     indexes (default model.rxns).
% useComps : logical, optional
%     include the compartment of each metabolite (default true).
% sortRevRxns : logical, optional
%     sort reversible reactions so that the metabolite that is first in the
%     lexicographic order is a reactant (default false).
% sortMetNames : logical, optional
%     sort the metabolite names in the equation. Uses compartment even if
%     useComps is false (default false).
% useMetID : logical, optional
%     use metabolite ID in generated equations (default false).
% useFormula : logical, optional
%     use metabolite formula in generated equations (default false).
% useRevField : logical, optional
%     use the model.rev field to indicate reaction reversibility,
%     alternatively this is determined from the model.ub and model.lb fields
%     (default true).
%
% Returns
% -------
% equationStrings : cell
%     a cell array with equations.
%
% Examples
% --------
%     equationStrings = constructEquations(model, rxns, useComps, ...
%         sortRevRxns, sortMetNames, useMetID, useFormula, useRevField);
%
% Notes
% -----
% If useRevField is false, then reactions should be organized in their
% forward direction (e.g. ub = 1000 and lb = -1000/0) for the reversibility
% to be correctly determined.

p=parseRAVENargs(varargin, {'rxns',[]; 'useComps',true; 'sortRevRxns',false; 'sortMetNames',false; 'useMetID',false; 'useFormula',false; 'useRevField',true});
rxns=p.rxns;
if isempty(rxns)
    rxns=model.rxns;
elseif ~islogical(rxns) && ~isnumeric(rxns)
    rxns=convertCharArray(rxns);
end
useComps=p.useComps;
sortRevRxns=p.sortRevRxns;
sortMetNames=p.sortMetNames;
useMetID=p.useMetID;
useFormula=p.useFormula;
useRevField=p.useRevField;

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
