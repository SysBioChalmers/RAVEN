function printModel(model,rxnList,outputString,outputFile,metaboliteList)
% printModel
%   Prints reactions to the screen or to a file
%
%   model           a model structure
%   rxnList         either a cell array of reaction IDs, a logical vector 
%                   with the same number of elements as reactions in the model,
%                   or a vector of indexes to remove (opt, default
%                   model.rxns)
%   outputString    a string that specifies the output of each reaction (opt,
%                   default '%rxnID (%rxnName)\n\t%eqn [%lower %upper]\n')
%   outputFile      a file to save the print-out to (opt, default is output to
%                   the command window)
%   metaboliteList  cell array of metabolite names. Only reactions
%                   involving any of these metabolites will be 
%                   printed (opt)
%
%   The following codes are available for user-defined output strings:
%
%   %rxnID      reaction ID
%   %rxnName    reaction name
%   %lower      lower bound
%   %upper      upper bound
%   %obj        objective coefficient
%   %eqn        equation
%   %element    equation using the metabolite formulas rather than
%               metabolite names
%   %unbalanced "(*)" if the reaction is unbalanced and "(-)" if it could not
%               be parsed
%   %lumped     equation where the elemental compositions for the left/right
%               hand sides are lumped
%
%   NOTE: This is just a wrapper function around printFluxes. It is
%           intended to be used when there is no flux distribution.
%
%   Usage: printModel(model,rxnList,outputString,outputFile,metaboliteList)

if nargin<2
    rxnList=model.rxns;
end
if isempty(rxnList)
    rxnList=model.rxns;
end
if nargin<3 || isempty(outputString)
    outputString='%rxnID (%rxnName)\n\t%eqn [%lower %upper]\n';
else
    outputString=char(outputString);
end
if nargin<4
    outputFile=[];
else
    outputFile=char(outputFile);
end
if nargin<5
    metaboliteList=[];
end

I=getIndexes(model,rxnList,'rxns',true)*1.00; %To convert it to "fluxes"

if ~isempty(metaboliteList)
    printFluxes(model, I, false, 0.1, outputFile,outputString,metaboliteList);
else
    printFluxes(model, I, false, 0.1, outputFile,outputString);
end
end
