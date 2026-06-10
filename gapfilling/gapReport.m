function [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat, canProduceWithoutInput, canConsumeWithoutOutput, ...
    connectedFromTemplates, addedFromTemplates]=gapReport(model, templateModels)
% gapReport  Perform a gap analysis and summarize the results.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% templateModels : cell, optional
%     a cell array of template models to use for gap filling.
%
% Returns
% -------
% noFluxRxns : cell
%     reactions that cannot carry flux.
% noFluxRxnsRelaxed : cell
%     reactions that cannot carry flux even if the mass balance
%     constraint is relaxed so that net production of all metabolites is
%     allowed.
% subGraphs : struct
%     the metabolites in each of the isolated sub-networks.
% notProducedMets : cell
%     the metabolites that could not have net production.
% minToConnect : struct
%     the minimal number of metabolites that need to be connected in
%     order to be able to produce all other metabolites, and which
%     metabolites each of them connects.
% neededForProductionMat : double
%     matrix where n x m is true if metabolite n allows for production of
%     metabolite m.
% canProduceWithoutInput : cell
%     metabolites that could be produced even when there is no input to
%     the model.
% canConsumeWithoutOutput : cell
%     metabolites that could be consumed even when there is no output from
%     the model.
% connectedFromTemplates : cell
%     the reactions that could be connected using the template models.
% addedFromTemplates : struct
%     the reactions that were added from the template models and which
%     model they were added from.
%
% Examples
% --------
%     [noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, ...
%         minToConnect, neededForProductionMat, connectedFromTemplates, ...
%         addedFromTemplates] = gapReport(model, templateModels);

if nargin<2
    templateModels=[];
    connectedFromTemplates=[];
    addedFromTemplates=[];
end

fprintf(['Gap analysis for ' model.id ' - ' model.name '\n\n']);
if isfield(model,'unconstrained')
    calculateINOUT=true;
    closedModel=model;
    model=simplifyModel(model);
else
    canConsumeWithoutOutput={};
    canProduceWithoutInput={};
    calculateINOUT=false;
end

model2=model;
model2.b=[model2.b inf(numel(model2.mets),1)];
I=haveFlux(model);
noFluxRxns=model.rxns(~I);
J=haveFlux(model2);
noFluxRxnsRelaxed=model2.rxns(~J);
bModel=removeReactions(model,~I,true,true);
cModel=removeReactions(model2,~J,true,true);
fprintf('***Overview\n');
fprintf([num2str(numel(model.rxns)-sum(I)) ' out of ' num2str(numel(model.rxns))...
    ' reactions cannot carry flux (' num2str(numel(model.rxns)-sum(J)) ' if net production of all metabolites is allowed)\n']);
fprintf([num2str(numel(model.mets)-numel(bModel.mets)) ' out of ' num2str(numel(model.mets))...
    ' metabolites are unreachable (' num2str(numel(model.mets)-numel(cModel.mets)) ' if net production of all metabolites is allowed)\n']);

fprintf('\n***Isolated subnetworks\n');
subGraphs=getAllSubGraphs(model);
fprintf(['A total of ' num2str(size(subGraphs,2)) ' isolated sub-networks are present in the model\n']);
for i=1:size(subGraphs,2)
    fprintf(['\t' num2str(i) '. ' num2str(sum(subGraphs(:,i))) ' metabolites\n']);
end

fprintf('\n***Metabolite connectivity\n');
[notProducedMets, ~, neededForProductionMat,minToConnect]=checkProduction(model,true,model.comps,false);
fprintf(['To enable net production of all metabolites, a total of ' num2str(numel(minToConnect)) ' metabolites must be connected\n']);
fprintf('Top 10 metabolites to connect:\n');
for i=1:min(10,numel(minToConnect))
    fprintf(['\t' num2str(i) '. ' minToConnect{i} '\n']);
end

if calculateINOUT==true
    fprintf('\n***Mass balancing\n');
    produced=canProduce(closedModel);
    canProduceWithoutInput=closedModel.mets(produced);
    consumed=canConsume(closedModel);
    canConsumeWithoutOutput=closedModel.mets(consumed);
    fprintf([num2str(numel(canConsumeWithoutOutput)) ' metabolites could be consumed without any outputs\n' num2str(numel(canProduceWithoutInput)) ' metabolites could be produced without any inputs\n']);
end

if ~isempty(templateModels)
    fprintf('\n***Automated gap-filling\n');
    [connectedFromTemplates, ~, addedFromTemplates]=fillGaps(model,templateModels);
    t=templateModels{1}.id;
    for i=2:numel(templateModels)
        t=[t ', ' templateModels{i}.id];
    end
    fprintf([num2str(numel(connectedFromTemplates)) ' unconnected reactions can be connected by including ' num2str(numel(addedFromTemplates)) ' reactions from\n' t '\n']);
end
end
