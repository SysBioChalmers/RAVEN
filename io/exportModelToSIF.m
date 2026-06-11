function exportModelToSIF(model,fileName,varargin)
% exportModelToSIF  Export a constraint-based model to a SIF file.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% fileName : char
%     the filename to export the model to.
% graphType : char, optional
%     the type of graph to export to (default 'rc'):
%
%     - 'rc' : reaction-compound
%     - 'rr' : reaction-reaction
%     - 'cc' : compound-compound
% rxnLabels : cell, optional
%     cell array with labels for reactions (default model.rxns).
% metLabels : cell, optional
%     cell array with labels for metabolites (default model.mets).
%
% Examples
% --------
%     exportModelToSIF(model, fileName, graphType, rxnLabels, metLabels);
fileName=char(fileName);
p=parseRAVENargs(varargin, {'graphType',[]; 'rxnLabels',[]; 'metLabels',[]});
graphType=p.graphType; rxnLabels=p.rxnLabels; metLabels=p.metLabels;
if isempty(graphType)
    graphType='rc';
else
    graphType=char(graphType);
end

if isempty(rxnLabels)
    rxnLabels=model.rxns;
else
    rxnLabels=convertCharArray(rxnLabels);
end
if isempty(metLabels)
    metLabels=model.mets;
else
    metLabels=convertCharArray(metLabels);
end

if ~strcmpi(graphType,'rc') && ~strcmpi(graphType,'rr') && ~strcmpi(graphType,'cc')
    EM='The graph type is incorrect';
    dispEM(EM);
end

if numel(rxnLabels)~=numel(unique(rxnLabels))
    EM='Not all reaction labels are unique';
    dispEM(EM,false);
end
if numel(metLabels)~=numel(unique(metLabels))
    EM='Not all metabolite labels are unique';
    dispEM(EM,false);
end

if strcmpi(graphType,'rc')
    G=model.S;
    A=rxnLabels;
    B=metLabels;
end
if strcmpi(graphType,'rr')
    G=model.S'*model.S;
    A=rxnLabels;
    B=rxnLabels;
end
if strcmpi(graphType,'cc')
    %A metabolite is linked to all products of the reactions that it
    %participates in If G=model.S*model.S' then all connections will be
    %double, which looks weird
    irrevModel=convertToIrrev(model);
    G=sparse(numel(model.mets),numel(model.mets));
    for i=1:numel(model.mets)
        I=irrevModel.S(i,:)<0; %Get the reactions in which it is a substrate
        [J, ~]=find(irrevModel.S(:,I)>0);
        G(J,i)=1;
    end
    A=metLabels;
    B=metLabels;
end

fid=fopen(fileName,'w');

for i=1:size(G,2)
    I=G(:,i)~=0;
    nodes=setdiff(B(I),A(i)); %Don't include connection to itself
    if ~isempty(nodes)
        nNodes=numel(nodes);
        nodes(1:nNodes-1)=strcat(nodes(1:nNodes-1),{'\t'});
        fullString=[nodes{:}];
        fprintf(fid,[A{i} '\t' graphType '\t' fullString '\n']);
    end
end
fclose(fid);
end
