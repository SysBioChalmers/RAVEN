function spontaneousRxnList=addSpontaneousRxns(rxnList, metList)
% addSpontaneousRxns
%   Retrieve spontaneous reactions based on the pathway-spontaneousRxn
%   associations curated by the MetaCyc database
%
%   rxnList              query list of reactions in cell array
%
%   metList              query list of metabolites in cell array
%
%   spontaneousRxnList   reterieved spontaneous reactions associated to
%                        the queried MetaCyc reactions and metabolites
%
%   Usage: spontaneousRxnList=addSpontaneousRxns(rxnList, metList)
%
%   Hao Wang, 2018-03-01
%

if nargin<2
    disp('Missing input arguments');
    return;
end

rxnList=unique(rxnList);

% Create the matrix of MetaCyc pathways and spontaneous reactions
load('metaCycRxns.mat');
pathways={};
for i=1:numel(metaCycRxns.rxns)
    if ~isempty(metaCycRxns.pwys{i})
        pathways=[pathways;transpose(strsplit(metaCycRxns.pwys{i},';'))];
    end
end
pathways=unique(pathways);

% Genearte the matirx, row: spontaneousRxns, column: pathways
sprxnPwyMat=zeros(numel(isSpontaneous),numel(pathways));
for i=1:numel(isSpontaneous)
    [a, b]=ismember(isSpontaneous{i},metaCycRxns.rxns);
    if ~isempty(metaCycRxns.pwys{b})
        [crap, indexes]=ismember(strsplit(metaCycRxns.pwys{b},';'),pathways);
        sprxnPwyMat(i,indexes)=1;   
	end
end

% Go through the rxnList and obtain relevant pathways
pwys={};
for i=1:numel(rxnList)
    [a, b]=ismember(rxnList{i},metaCycRxns.rxns);
	if a && ~isempty(metaCycRxns.pwys{b})
        pwys=[pwys;transpose(strsplit(metaCycRxns.pwys{b},';'))];
    end
end
pwys=unique(pwys);

% Get spontaneous reactions associated with the pathways
hits=[];
for i=1:numel(pwys)
    [a, b]=ismember(pwys{i},pathways);
    if ~isempty(find(sprxnPwyMat(:,b)))
        hits=[hits;find(sprxnPwyMat(:,b))];
    end		
end
spontaneousRxnList=isSpontaneous(unique(hits));

% Check if the reactants/products are included in metList
for i=1:numel(spontaneousRxnList)
    [a, b]=ismember(spontaneousRxnList{i},metaCycRxns.rxns);

    %Obtain the reactants and products
    reactants=all(ismember(metaCycRxns.mets(find(metaCycRxns.S(:,b)==-1)),metList));
    products=all(ismember(metaCycRxns.mets(find(metaCycRxns.S(:,b)==1)),metList));

    %Either reactants or products should be present for reversible reaction
    if metaCycRxns.rev(b)==1
        if ~(reactants || products)
            spontaneousRxnList{i}=[];
        end
    %Reactants should be present for irreversible reactions
    else
        if ~reactants
            spontaneousRxnList{i}=[];
        end	
    end
end
spontaneousRxnList=spontaneousRxnList(~cellfun(@isempty, spontaneousRxnList));

% Remove the spontaneousRxns that have already been included
spontaneousRxnList=setdiff(spontaneousRxnList,rxnList);

end
