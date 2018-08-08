function compStruct=compareModels(models,printResults)
% compareModels
%   Compares two or more models with respect to overlap in terms of genes,
%   reactions, metabolites and compartments.
%
%   models              cell array of two or more models
%   printResults        true if the results should be printed on the screen
%                       (opt, default false)
%
%   compStruct          structure that contains the comparison
%       modelIDs        cell array of model ids
%       rxns            These contain the comparison for each field. 'equ' are
%                       the equations after sorting and 'uEqu' are the
%                       equations when not taking compartmentalization into acount 
%       mets
%       genes
%       eccodes
%       metNames
%       equ
%       uEqu
%           comparison	binary matrix where each row indicate which models are
%                       included in the comparison
%           nElements   vector with the number of elements for each
%                       comparison
%
%   Usage: compStruct=compareModels(models,printResults)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    printResults=true;
end

if numel(models)<=1
   dispEM('Cannot compare only one model. Use printModelStats if you want a summary of a model'); 
end

compStruct.modelIDs={};
for i=1:numel(models)
    compStruct.modelIDs=[compStruct.modelIDs;models{i}.id];
    models{i}.equ=constructEquations(models{i},models{i}.rxns,true,true,true);
    models{i}.uEqu=constructEquations(models{i},models{i}.rxns,false,true,true);
end

field='rxns';
compStruct.rxns.comparison=getToCheck(models,field);
compStruct.rxns.nElements=checkStuff(getElements(models,field),compStruct.rxns.comparison);
if printResults==true
    fprintf('*** Comparison of reaction IDs:\n');
    printList(models,compStruct.rxns.comparison,compStruct.rxns.nElements);
    fprintf('\n\n');
end

field='mets';
compStruct.mets.comparison=getToCheck(models,field);
compStruct.mets.nElements=checkStuff(getElements(models,field),compStruct.mets.comparison);
if printResults==true
    fprintf('*** Comparison of metabolite IDs:\n');
    printList(models,compStruct.mets.comparison,compStruct.mets.nElements);
    fprintf('\n\n');
end

field='genes';
compStruct.genes.comparison=getToCheck(models,field);
compStruct.genes.nElements=checkStuff(getElements(models,field),compStruct.genes.comparison);
if printResults==true
    fprintf('*** Comparison of gene IDs:\n');
    printList(models,compStruct.genes.comparison,compStruct.genes.nElements);
    fprintf('\n\n');
end

field='eccodes';
compStruct.eccodes.comparison=getToCheck(models,field);
compStruct.eccodes.nElements=checkStuff(getElements(models,field),compStruct.eccodes.comparison);
if printResults==true
    fprintf('*** Comparison of ec-numbers:\n');
    printList(models,compStruct.eccodes.comparison,compStruct.eccodes.nElements);
    fprintf('\n\n');
end

field='metNames';
compStruct.metNames.comparison=getToCheck(models,field);
compStruct.metNames.nElements=checkStuff(getElements(models,field),compStruct.metNames.comparison);
if printResults==true
    fprintf('*** Comparison of metabolite names:\n');
    printList(models,compStruct.metNames.comparison,compStruct.metNames.nElements);
    fprintf('\n\n');
end

field='equ';
compStruct.equ.comparison=getToCheck(models,field);
compStruct.equ.nElements=checkStuff(getElements(models,field),compStruct.equ.comparison);
if printResults==true
    fprintf('*** Comparison of equations with compartment:\n');
    printList(models,compStruct.equ.comparison,compStruct.equ.nElements);
    fprintf('\n\n');
end

field='uEqu';
compStruct.uEqu.comparison=getToCheck(models,field);
compStruct.uEqu.nElements=checkStuff(getElements(models,field),compStruct.uEqu.comparison);
if printResults==true
    fprintf('*** Comparison of equations without compartment:\n');
    printList(models,compStruct.uEqu.comparison,compStruct.uEqu.nElements);
    fprintf('\n\n');
end
end
function A=getElements(models,field)
    A={};
    for i=1:numel(models)
       if isfield(models{i},field) 
        A=[A;{getfield(models{i},field)}]; 
       end
    end
end
function toCheck=getToCheck(models,field)
    %Get all the combinations that should be checked for overlap (including the
    %single ones)
    toCheckA=[];
    I=find(cellfun(@checkField,models));
    nI=numel(I);
    for i=nI:-1:1
        combs=combnk(1:nI,i);
        toAdd=false(size(combs,1),nI);
        for j=1:size(combs,1)
           toAdd(j,combs(j,:))=true; 
        end
        toCheckA=[toCheckA;toAdd];
    end
    
    %If not all of the models have the required field 
    toCheck=false(size(toCheckA,1),numel(models));
    toCheck(:,I)=toCheckA;
    
    %Ugly thing to get around parameters
    function I=checkField(A)
        I=isfield(A,field);
    end
end
function printList(models,toCheck,nElements)
    %To guess how many spaces that are needed to align
    firstLen=[];
    for i=1:size(toCheck,1)
       label=[];
       I=find(toCheck(i,:));
       for j=1:numel(I)
            label=[label models{I(j)}.id '/'];
       end
       if i==1
           firstLen=numel(label);
       end
       nSpaces=firstLen-numel(label);
       fprintf([label(1:end-1) '  ' repmat(sprintf(' '),1,nSpaces) num2str(nElements(i)) '\n']);
    end
end
function nElements=checkStuff(A,toCheck)
    %Now loop through the toCheck matrix, starting with the combination with the
    %most models. Only elements that weren't in iteration n are considered in
    %iteration n+1.
    nElements=zeros(size(toCheck,1),1);
    alreadyChecked=[];
    for i=1:size(toCheck,1)
        I=find(toCheck(i,:));
        inCommon=setdiff(A{I(1)},alreadyChecked);
        for j=2:numel(I)
           inCommon=intersect(inCommon,A{I(j)}); 
        end
        alreadyChecked=union(alreadyChecked,inCommon);
        nElements(i)=numel(inCommon);
    end
end