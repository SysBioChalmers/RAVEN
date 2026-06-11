function printModel(model,varargin)
% printModel  Print reactions to the screen or to a file.
%
% This is a wrapper around printFluxes, intended for use when there is no
% flux distribution.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% rxnList : cell or logical or double
%     either a cell array of reaction IDs, a logical vector with the same
%     number of elements as reactions in the model, or a vector of indexes
%     to print (default model.rxns).
% outputString : char
%     a string that specifies the output of each reaction (default
%     '%rxnID (%rxnName)\n\t%eqn [%lower %upper]\n').
% outputFile : char
%     a file to save the print-out to (default is output to the command
%     window).
% metaboliteList : cell
%     cell array of metabolite names. Only reactions involving any of these
%     metabolites will be printed.
%
% Notes
% -----
% The following codes are available for user-defined output strings:
%
% - %rxnID : reaction ID
% - %rxnName : reaction name
% - %lower : lower bound
% - %upper : upper bound
% - %obj : objective coefficient
% - %eqn : equation
% - %element : equation using the metabolite formulas rather than metabolite
%   names
% - %unbalanced : "(*)" if the reaction is unbalanced and "(-)" if it could
%   not be parsed
% - %lumped : equation where the elemental compositions for the left/right
%   hand sides are lumped
%
% Examples
% --------
%     printModel(model, rxnList, outputString, outputFile, metaboliteList);

p=parseRAVENargs(varargin, {'rxnList',[]; 'outputString',[]; 'outputFile',[]; 'metaboliteList',[]});
rxnList=p.rxnList;
if isempty(rxnList)
    rxnList=model.rxns;
elseif ~islogical(rxnList) && ~isnumeric(rxnList)
    rxnList=convertCharArray(rxnList);
end
outputString=p.outputString;
if isempty(outputString)
    outputString='%rxnID (%rxnName)\n\t%eqn [%lower %upper]\n';
else
    outputString=char(outputString);
end
outputFile=p.outputFile;
if ~isempty(outputFile)
    outputFile=char(outputFile);
end
metaboliteList=p.metaboliteList;
if ~isempty(metaboliteList)
    metaboliteList=convertCharArray(metaboliteList);
end

I=getIndexes(model,rxnList,'rxns',true)*1.00; %To convert it to "fluxes"

if ~isempty(metaboliteList)
    printFluxes(model, I, false, 0.1, outputFile,outputString,metaboliteList);
else
    printFluxes(model, I, false, 0.1, outputFile,outputString);
end
end
