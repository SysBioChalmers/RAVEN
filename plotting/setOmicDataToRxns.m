function [v1] = setOmicDataToRxns(omics, model)
% USAGE:
% [v1] = setOmicDataToRxns(transcriptome, model)
% INPUTS:
% transcriptome         A two columns entrace with gene names and
%                       log-fold expression changes
% model                 A model with a COBRA structure
% OUTPUTS:
% v1                    A one column entrace with expression changes in
%                       metabolic genes.
    j = 1;
    for i = 1:length(model.genes)
        a = find(ismember(omics(:,1), model.genes{i}));
        if ~isempty(a)
            genesList{j,1} = omics{a,1};
            genesList{j,2} = omics{a,2};
            j = j + 1;
        end
    end
    rxnsFromGenes = findRxnsFromGenes(model, genesList(:,1));
    s = struct2cell(rxnsFromGenes);
    k = 1;
    for i = 1:length(s)
        s1 = s{i};
        [l,~]=size(s1);
        for j = 1:l
            s2=s1{j};
            rxnsList{k,1}=s2;
            k = k + 1;
        end
    end
    rxnsList = (unique(rxnsList));
    
    for i = 1:length(rxnsList)
        rxnsList {i,2} = findGenesFromRxns(model, rxnsList{i});
    end
    v1 = zeros(length(model.rxns),1);
    for i = 1:length(rxnsList)
        rxnGenes = rxnsList{i,2};
        indRxnGene = rxnGenes{1,1};
        c = 0.0;
        for j = 1:length(indRxnGene)
            a = find(ismember(genesList(:,1),indRxnGene{j,1}));
            if ~isempty(a) && abs(cell2mat(genesList(a,2))) >= c
                c = cell2mat(genesList(a,2));
            end
        end
        rxnID = findRxnIDs(model,rxnsList{i,1});
        v1(rxnID,1) = c;
    end
end