function printModelStats(model, printModelIssues, printDetails)
% printModelStats
%   prints some statistics about a model to the screen
%
%   model               a model structure
%   printModelIssues    true if information about unconnected
%                       reactions/metabolites and elemental balancing
%                       should be printed (opt, default false)
%   printDetails        true if detailed information should be printed
%                       about model issues. Only used if printModelIssues
%                       is true (opt, default true)
%
%   Usage: printModelStats(model,printModelIssues, printDetails)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    printModelIssues=false;
end
if nargin<3
    printDetails=true;
end

fprintf(['Network statistics for ' model.id ': ' model.description '\n']);

%Get which reactions are present in each compartment
rxnComps=sparse(numel(model.rxns),numel(model.comps));

%For each compartment, find the metabolites that are present in that
%compartment and then the reactions they are involved in
for i=1:numel(model.comps)
    [~, I]=find(model.S(model.metComps==i,:));
    rxnComps(I,i)=1;
end

if isfield(model,'eccodes')
    fprintf(['EC-numbers\t\t\t' num2str(numel(unique(model.eccodes))) '\n']);
end

%Print information about genes
if isfield(model,'genes')
    fprintf(['Genes*\t\t\t\t' num2str(numel(model.genes)) '\n']);
    %Find the genes in each compartment
    for i=1:numel(model.comps)
        [~, I]=find(model.rxnGeneMat(rxnComps(:,i)==1,:));
        fprintf(['\t' model.compNames{i} '\t' num2str(numel(unique(I))) '\n']);
    end
end

%Print information about reactions
fprintf(['\nReactions*\t\t\t' num2str(numel(model.rxns)) '\n']);
for i=1:numel(model.comps)
    fprintf(['\t' model.compNames{i} '\t' num2str(sum(rxnComps(:,i))) '\n']);
end

%Removes the effect of compartments and removes duplicate reactions
temp=model;
temp.comps(:)={'s'}; %Set all compartments to be the same
equ=constructEquations(sortModel(temp,true,true),temp.rxns,false);

fprintf(['Unique reactions**\t' num2str(numel(unique(equ))) '\n']);

%Print information about metabolites
fprintf(['\nMetabolites\t\t\t' num2str(numel(model.mets)) '\n']);
for i=1:numel(model.comps)
    fprintf(['\t' model.compNames{i} '\t' num2str(sum(model.metComps==i)) '\n']);
end

fprintf(['Unique metabolites\t' num2str(numel(unique(model.metNames))) '\n']);

fprintf('\n* Genes and reactions are counted for each compartment if any of the corresponding metabolites are in that compartment. The sum may therefore not add up to the total number.\n');
fprintf('** Unique reactions are defined as being biochemically unique (no compartmentalization)\n');

%Also print some potential problems if there are any
if printModelIssues==true
    fprintf(['\nShort model quality summary for ' model.id ': ' model.description '\n']);
    
    %Check that all the metabolites are being used
    involvedMat=model.S;
    involvedMat(involvedMat~=0)=1;
    usedMets=sum(involvedMat,2);
    notPresent=find(usedMets==0);
    if ~isempty(notPresent)
        errorText=['Non-used metabolites\t' num2str(numel(notPresent)) '\n'];
        if printDetails==true
            for i=1:numel(notPresent)
                errorText=[errorText '\t(' model.mets{notPresent(i)} ') ' model.metNames{notPresent(i)} '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
    
    %Check if there are empty reactions
    usedRxns=sum(involvedMat,1);
    notUsed=find(usedRxns==0);
    if ~isempty(notUsed)
        errorText=['Empty reactions\t' num2str(numel(notUsed)) '\n'];
        if printDetails==true
            for i=1:numel(notUsed)
                errorText=[errorText '\t' model.rxns{notUsed(i)} '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
    
    %Check if there are dead-end reactions/metabolites
    [~, deletedReactions, deletedMetabolites]=simplifyModel(model,true,false,false,true);
    
    if ~isempty(deletedReactions)
        errorText=['Dead-end reactions\t' num2str(numel(deletedReactions)) '\n'];
        if printDetails==true
            for i=1:numel(deletedReactions)
                errorText=[errorText '\t' deletedReactions{i} '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
    
    %Ignore non-used metabolites
    deletedMetabolites=setdiff(deletedMetabolites,model.mets(notPresent));
    %Must map to indexes in order to print names
    deletedMetabolites=find(ismember(model.mets,deletedMetabolites));
    if ~isempty(deletedMetabolites)
        errorText=['Dead-end metabolites\t' num2str(numel(deletedMetabolites)) '\n'];
        if printDetails==true
            for i=1:numel(deletedMetabolites)
                errorText=[errorText '\t(' model.mets{deletedMetabolites(i)} ') ' model.metNames{deletedMetabolites(i)} '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
    
    balanceStructure=getElementalBalance(model);
    
    notParsed=find(balanceStructure.balanceStatus<0);
    notBalanced=find(balanceStructure.balanceStatus==0);
    
    if ~isempty(notParsed)
        errorText=['Reactions which could not be elementally balanced\t' num2str(numel(notParsed)) '\n'];
        if printDetails==true
            for i=1:numel(notParsed)
                errorText=[errorText '\t' model.rxns{notParsed(i)} '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
    if ~isempty(notBalanced)
        errorText=['Reactions which are elementally unbalanced\t' num2str(numel(notBalanced)) '\n'];
        if printDetails==true
            names=strcat(balanceStructure.elements.names,{', '});
            for i=1:numel(notBalanced)
                badOnes=sprintf('%s', names{abs(balanceStructure.leftComp(notBalanced(i),:)-balanceStructure.rightComp(notBalanced(i),:))>10^-7});
                errorText=[errorText '\t' model.rxns{notBalanced(i)} '\t' badOnes(1:end-2) '\n'];
            end
            errorText=[errorText '\n'];
        end
        fprintf(errorText);
    end
end
end
