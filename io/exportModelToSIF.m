function exportModelToSIF(model,fileName,graphType,rxnLabels,metLabels)
% exportModelToSIF
%   Exports a constraint-based model to a SIF file
%
%   model         a model structure
%   fileName      the filename  to export the model to
%   graphType     the type of graph to export to (opt, default 'rc')
%                 'rc'  reaction-compound
%                 'rr'  reaction-reaction
%                 'cc'  compound-compound
%   rxnLabels     cell array with labels for reactions (opt, default
%                 model.rxns)
%   metLabels     cell array with labels for metabolites (opt, default
%                 model.mets)
%
%   Usage: exportModelToSIF(model,fileName,graphType,rxnLabels,metLabels)
%
%   Rasmus Agren, 2014-01-08
%

if nargin<3
   graphType='rc';
end

if nargin<4
   rxnLabels=model.rxns;
end
if nargin<5
   metLabels=model.mets;
end
if isempty(rxnLabels)
   rxnLabels=model.rxns;
end
if isempty(metLabels)
   metLabels=model.mets;
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
   %A metabolite is linked to all products of the reactions that it participates in
   %If G=model.S*model.S' then all connections will be double, which looks
   %weird
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
