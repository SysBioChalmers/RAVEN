function equationStrings=constructEquations(model,rxns,useComps,sortRevRxns,sortMetNames)
% constructEquations
%   Construct equation strings for reactions
%
%   model             a model structure
%   rxns              either a cell array of reaction IDs, a logical vector with the 
%                     same number of elements as reactions in the model, or a vector 
%                     of reaction indexes (opt, default model.rxns)
%	useComps          include the compartment of each metabolite (opt, default true)
%	sortRevRxns       sort reversible reactions so that the metabolite that is first in
%                     the lexiographic order is a reactant (opt, default
%                     false)
%	sortMetNames      sort the metabolite names in the equation. Uses
%                     compartment even if useComps is false (opt, default
%                     false)
%
%   equationStrings     a cell array with equations
%
% 	Usage: equationStrings=constructEquations(model,rxns,useComps,...
%           sortRevRxns,sortMetNames)
%
%   Rasmus Agren, 2013-05-16
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
if isempty(rxns)
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
    reactants=find(model.S(:,indexes(i))<0);
    products=find(model.S(:,indexes(i))>0);
    eqn='';

    for j=1:numel(reactants)
    	if j==1
        	plusString='';
        else
        	plusString=' + ';
        end

        stoich=num2str(model.S(reactants(j),indexes(i))*-1);
        
        if str2double(stoich)==1
        	stoich='';
        else
        	stoich=[stoich ' '];
        end
        
        if useComps==true
            eqn=[eqn plusString stoich model.metNames{reactants(j)} '[' model.comps{model.metComps(reactants(j))} ']'];
        else
            eqn=[eqn plusString stoich model.metNames{reactants(j)}];
        end
	end
       
    if model.rev(indexes(i))==0
    	eqn=[eqn ' => '];
    else
    	eqn=[eqn ' <=> '];
    end
       
    for j=1:numel(products)
    	if j==1
        	plusString='';
        else
        	plusString=' + ';
        end

        stoich=num2str(model.S(products(j),indexes(i)));
        
        if str2double(stoich)==1
        	stoich='';
        else
        	stoich=[stoich ' '];
        end
        
        if useComps==true
            eqn=[eqn plusString stoich model.metNames{products(j)} '[' model.comps{model.metComps(products(j))} ']'];
        else
            eqn=[eqn plusString stoich model.metNames{products(j)}];
        end
    end
    equationStrings{i}=eqn;
end
end
