function printFluxes(model, fluxes, varargin)
% printFluxes  Print reactions and fluxes to the screen or to a file.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% fluxes : double
%     a vector with fluxes.
% onlyExchange : logical, optional
%     only print exchange fluxes (default true).
% cutOffFlux : double, optional
%     only print fluxes with absolute values above or equal to this value
%     (default 10^-8).
% outputFile : char, optional
%     a file to save the print-out to (default is output to the command
%     window).
% outputString : char, optional
%     a string that specifies the output of each reaction (default
%     '%rxnID\t(%rxnName):\t%flux\n').
% metaboliteList : cell, optional
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
% - %flux : flux
% - %element : equation using the metabolite formulas rather than metabolite
%   names
% - %unbalanced : "(*)" if the reaction is unbalanced and "(-)" if it could
%   not be parsed
% - %lumped : equation where the elemental compositions for the left/right
%   hand sides are lumped
%
% Examples
% --------
%     printFluxes(model, fluxes, onlyExchange, cutOffFlux, outputFile, ...
%                 outputString, metaboliteList);

p=parseRAVENargs(varargin, {'onlyExchange',true; 'cutOffFlux',10^-8; 'outputFile',[]; 'outputString',[]; 'metaboliteList',{}});
onlyExchange=p.onlyExchange;
cutOffFlux=p.cutOffFlux;
if isempty(cutOffFlux)
    cutOffFlux=10^-8;
end
outputFile=p.outputFile;
if ~isempty(outputFile)
    outputFile=char(outputFile);
    fid=fopen(outputFile,'w');
else
    fid=1;
end
outputString=p.outputString;
if isempty(outputString)
    outputString='%rxnID\t(%rxnName):\t%flux\n';
else
    outputString=char(outputString);
end
metaboliteList=p.metaboliteList;
if ~isempty(metaboliteList)
    metaboliteList=convertCharArray(metaboliteList);
end
if isempty(fluxes)
    EM='Empty vector of fluxes, solveLP possibly returned infeasible';
    dispEM(EM);
elseif size(fluxes,1)~=numel(model.rxns)
    EM='The number of fluxes and the number of reactions must be the same';
    dispEM(EM);
end

%Only keep reactions involving the defined metabolites
if ~isempty(metaboliteList)
    I=ismember(upper(model.metNames),upper(metaboliteList));
    [~, K]=find(model.S(I,:));
    
    %Delete all other reactions
    toDelete=true(numel(model.rxns),1);
    toDelete(K)=false;
    model=removeReactions(model,toDelete);
    fluxes(toDelete,:)=[];
end

if onlyExchange==true
    fprintf(fid,'EXCHANGE FLUXES:\n');
else
    fprintf(fid,'FLUXES:\n');
end

%Remove reactions which are below the cut off
toDelete=abs(fluxes)<cutOffFlux;
toDelete=all(toDelete,2);
model=removeReactions(model,toDelete,true,true);
fluxes(toDelete,:)=[];

if any(strfind(outputString,'%eqn'))
    %Construct the equations
    eqn=constructEquations(model);
else
    eqn=cell(numel(model.rxns),1);
    eqn(:)={''};
end
if any(strfind(outputString,'%element'))
    %For printing equations using the composition
    cModel=model;
    cModel.metNames=cModel.metFormulas;
    cModel.metNames(cellfun(@isempty,cModel.metNames))={'?'};
    element=constructEquations(cModel);
else
    element=cell(numel(model.rxns),1);
    element(:)={''};
end

if any(strfind(outputString,'%unbalanced')) || any(strfind(outputString,'%lumped'))
    balanceStructure=getElementalBalance(model);
end

unbalanced=cell(numel(model.rxns),1);
unbalanced(:)={''};
if any(strfind(outputString,'%unbalanced'))
    unbalanced(balanceStructure.balanceStatus==0)={'(*)'};
    unbalanced(balanceStructure.balanceStatus<0)={'(-)'};
end

lumped=cell(numel(model.rxns),1);
lumped(:)={''};
if any(strfind(outputString,'%lumped'))
    for i=1:numel(model.rxns)
        leftGroup='';
        rightGroup='';
        for j=1:numel(balanceStructure.elements.names)
            I=balanceStructure.leftComp(i,j);
            if I~=0
                if I==1
                    leftGroup=[leftGroup balanceStructure.elements.abbrevs{j}];
                else
                    leftGroup=[leftGroup balanceStructure.elements.abbrevs{j} num2str(I)];
                end
            end
            I=balanceStructure.rightComp(i,j);
            if I~=0
                if I==1
                    rightGroup=[rightGroup balanceStructure.elements.abbrevs{j}];
                else
                    rightGroup=[rightGroup balanceStructure.elements.abbrevs{j} num2str(I)];
                end
            end
        end
        if model.rev(i)
            lumped{i}=[leftGroup ' <=> ' rightGroup];
        else
            lumped{i}=[leftGroup ' => ' rightGroup];
        end
    end
end

for i=1:numel(model.rxns)
    %Only print if it's an exchange reaction or if all reactions should be
    %printed. Exchange reactions only have reactants or only products.
    reactants=model.S(:,i)<0;
    products=model.S(:,i)>0;
    
    %Only print if the absolute value is >= cutOffFlux
    if (onlyExchange==false || (~any(reactants) || ~any(products)))
        printString=outputString;
        
        %Produce the final string
        printString=strrep(printString,'%rxnID',model.rxns{i});
        printString=strrep(printString,'%eqn',eqn{i});
        printString=strrep(printString,'%rxnName',model.rxnNames{i});
        printString=strrep(printString,'%lower',num2str(model.lb(i)));
        printString=strrep(printString,'%upper',num2str(model.ub(i)));
        printString=strrep(printString,'%obj',num2str(model.c(i)));
        printString=strrep(printString,'%flux',num2str(fluxes(i,:)));
        printString=strrep(printString,'%element',element{i});
        printString=strrep(printString,'%unbalanced',unbalanced{i});
        printString=strrep(printString,'%lumped',lumped{i});
        fprintf(fid,printString);
    end
end

if fid~=1
    fprintf('File successfully saved.\n');
    fclose(fid);
end
end
