function spontaneousRxnList=retrieveSpontaneous(rxnList,metList)
% retrieveSpontaneous
%   Retrieve spontaneous reactions associated with the quried reactions
%
%   model                a model structure (will be the only argument)
%
%   rxnList              query list of reactions in cell array
%
%   metList              query list of metabolites in cell array
%
%   spontaneousRxnList   reterieved spontaneous reactions associated to
%                        the model (represented by rxnList and metList)
%
%   Usage: spontaneousRxnList=retrieveSpontaneous(rxnList,metList)
%
%   Hao Wang, 2017-07-23
%

if nargin<2
    disp('Missing input arguments');
    return;
end

rxnList=unique(rxnList);

% Create the matrix of spontaneous reaction and pathways
load('metaCycRxns.mat');
pathways={};
for i=1:numel(metaCycRxns.rxns)
		if ~isempty(metaCycRxns.pwys{i})
				pathways=[pathways;transpose(strsplit(metaCycRxns.pwys{i},';'))];
		end
end
pathways=unique(pathways);
sprxnPwyMat=zeros(numel(isSpontaneous),numel(pathways)); % row: sprxn, column: pathway
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
		if metaCycRxns.rev(b)==1
				if ~(reactants || products)
						spontaneousRxnList{i}=[];
				end
		else
				if ~reactants
						spontaneousRxnList{i}=[];
				end	
		end
end
spontaneousRxnList=spontaneousRxnList(~cellfun(@isempty, spontaneousRxnList));

% Remove the ones already included in the model
spontaneousRxnList=setdiff(spontaneousRxnList,rxnList);

end
