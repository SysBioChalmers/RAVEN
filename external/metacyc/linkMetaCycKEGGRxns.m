function linkMetaCycKEGGRxns
% linkMetaCycKEGGRxns
%   Link additional MetaCyc and KEGG reactions through metabolite mapping
%   This function only need to run once when the MetaCyc database is updated
%
%   NOTE: No arguments are required
%
%   Usage: linkMetaCycKEGGRxns

load('metaCycRxns.mat'); %load MetaCyc reactions
fprintf('NOTE: Importing MetaCyc reactions...\n');
metaCycModel = metaCycRxns;
metaCycModel.rxnFrom=cell(numel(metaCycModel.rxns),1);
metaCycModel.rxnFrom(:)={'MetaCyc'};
metaCycModel.grRules={};
metaCycModel.genes={};
metaCycModel.rxnGeneMat=sparse(numel(metaCycModel.rxns),1);


keggModel=getRxnsFromKEGG(); %load KEGG reactions

%Remove KEGG reactions with MetaCyc links based on original rxnLinks
rxnToRemove=intersect(rxnLinks.kegg,keggModel.rxns);
rxnToRemove=unique(rxnToRemove);
shrinkedKeggModel=removeReactions(keggModel,rxnToRemove,true,true);

fprintf('Mapping MetaCyc and KEGG reactions...\n');
%Resolve the shared but unmapped reactions (through mapping the involved
%metabolites). Replace mets information in KEGG model with the
%corresponding ones in MetaCyc.
load('metaCycMets.mat');
for i=1:numel(shrinkedKeggModel.mets)
    [a, b]=ismember(shrinkedKeggModel.mets{i},metaCycMets.keggid);
    if a
        %Replace the met ids
        shrinkedKeggModel.mets{i}=metaCycMets.mets{b};
        shrinkedKeggModel.metNames{i}=metaCycMets.mets{b};
    end
end

%Prepare for the merge of KEGG and MetaCyc super models

%Adding fields (comps, compNames, metNames, metComps)
metaCycModel.comps={'s'};
metaCycModel.compNames={'System'};
metaCycModel.metNames=metaCycModel.mets;
if ~isfield(metaCycModel,'metComps')
    metaCycModel.metComps=ones(numel(metaCycModel.mets),1);
end

shrinkedKeggModel.comps={'s'};
shrinkedKeggModel.compNames={'System'};
shrinkedKeggModel.metNames=shrinkedKeggModel.mets;
if ~isfield(shrinkedKeggModel,'metComps')
    shrinkedKeggModel.metComps=ones(numel(shrinkedKeggModel.mets),1);
end

%Merge models
mappingModel=mergeModels({shrinkedKeggModel metaCycModel},'metNames');

%Remove compounds proton and water because KEGG reactions often miss them
mappingModel=removeMets(mappingModel,{'PROTON','WATER'});

%Find same/similiar reactions by check identical columns in S matrix Go
%through the pure KEGG model and collect rxns id that have identical
%reactions to the ones in MetaCyc
I=find(ismember(mappingModel.rxnFrom,'KEGG'));
for i=1:numel(I)
    testRxn=mappingModel.S(:,I(i));
    match=find(~any(bsxfun(@minus, mappingModel.S,testRxn))); %Figure out the matched columns in S matrix
    % if ~isequal(match,I(i))
    count=numel(match);
    if count>1 %Only consider one-to-one match here
        for j=2:count
            %Only consider the match between KEGG and MetaCyc
            if isequal(mappingModel.rxnFrom{match(j)},'MetaCyc')
                %disp([mappingModel.rxns(match(1))
                %mappingModel.rxns(match(j))]); %keep this for later
                %checking
                rxnLinks.kegg=[rxnLinks.kegg;mappingModel.rxns(match(1))];
                rxnLinks.metacyc=[rxnLinks.metacyc;mappingModel.rxns(match(j))];
            end
        end
    end
    % end
end

% Load orginal rxnLinks
numLink=numel(rxnLinks.kegg);
rxnLinks.check=cell(numLink,1);
for i=1:numLink
    rxnLinks.check{i}=strcat(rxnLinks.kegg{i},rxnLinks.metacyc{i});
end

[~, index]=unique(rxnLinks.check);
rxnLinks.kegg=rxnLinks.kegg(index);
rxnLinks.metacyc=rxnLinks.metacyc(index);
rxnLinks=rmfield(rxnLinks,'check');

%Get the MetaCyc path and update the metaCycRxns.mat
ravenPath=findRAVENroot();
rxnsFile=fullfile(ravenPath,'external','metacyc','metaCycRxns.mat');
save(rxnsFile,'metaCycRxns','rxnLinks','TRANSPORT','UNBALANCED','UNDETERMINED','isSpontaneous');
fprintf(['Reaction associations between MetaCyc and KEGG have been successfully updated!\n\n']);

end
