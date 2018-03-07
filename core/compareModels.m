function compStruct=compareModels(models,printResults,plotResults,groupVector,fluxMets,fluxValues,objectiveFunctions,numRuns)
% compareModels
%   Compares two or more models with respect to overlap in terms of genes,
%   reactions, metabolites and compartments. New functionality also
%   compares model structural and functional similarity using
%   high-dimensional comparisons.
%
%   models              cell array of two or more models
%   printResults        true if the results should be printed on the screen
%                       (opt, default false)
%   plotResults         true if the results should be plotted
%                       (opt, default false)
%   groupVector         numeric vector for grouping similar models, i.e. by
%                       tissue (opt, default, all models ungrouped)
%   functionalComp      true if functional model comparison should be run
%                       (opt, default, false)
%   fluxMets            string array of metabolites used to constrain model
%                       simulations (opt, default = no constraints)
%   fluxValues          numeric vector (type: double) of values for
%                       constraining fluxMets during model simulations
%                       (opt, default = -1000, 1000)
%   objectiveFunctions  string array containing the objective function name
%                       for each model (opt, default is use objective
%                       function pre-set in model)
%   numRuns             single value selecting the number of simulations to
%                       perform for functional similarity comparisons
%                       (default, 30)
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
%       subsystems    
%           matrix      table with comparison of number of rxns per
%                       subsystem
%           ID          vector consisting of names of all subsystems
%       structComp      matrix with pairwise comparisons of model structure
%                       based on (1-Hamming distance) between models
%       structCompMap   matrix with 3D tSNE mapping of model structures
%                       based on Hamming distances
%       funcCompMap     matrix with 3D tSNE mapping of model fluxes using numRuns
%                       perturbations to the fluxMet input values "fluxValues" 
%                       (+/- 10% nominal) based on Euclidean distances
%
%   Usage: compStruct=compareModels(models,printResults,plotResults,...
%                     groupVector,fluxMets,fluxValues,objectiveFunctions,numRuns)
%
%   Rasmus Agren, 2014-02-07
%
%   Edited: Daniel Cook, 2018-03-12

%% Set up input defaults
if nargin<2
    printResults=false;
end

if nargin<3
    plotResults=false;
end

if nargin<4
    groupVector = [];
end

if nargin<5
    functionalComp = false;
end

if nargin<6
    fluxMets = [];
end

if nargin<7
    fluxValues = [];
end

if nargin<8
    objectiveFunction = [];
end

if nargin<9
    numRuns = 30;
end

if numel(models)<=1
    EM='Cannot compare only one model. Use printModelStats if you want a summary of a model';
    dispEM(EM);
end

%% Compare models (original)
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

%% Compare models structure & function based on high-dimensional methods
field = 'subSystems';
compStruct.subsystems.ID = catModelElements(models,field);
compStruct.subsystems.matrix = compareSubsystems(models,field);
if printResults==true
    % This could use come cleaning up
    fprintf('*** Comparison of reaction subsystem populations:\n');    
    str = [" "];
    for i=1:length(compStruct.modelIDs)
        str(1,i+1) = compStruct.modelIDs{i};
    end
    for i = 1:size(compStruct.subsystems.matrix,1)
        if length(compStruct.subsystems.ID{i}) <= 16
            str(i+1,1) = compStruct.subsystems.ID{i};
        else
            str(i+1,1) = [compStruct.subsystems.ID{i}(1:16) ' ...'];
        end
        for j=1:length(compStruct.modelIDs)
            str(i+1,j+1) = num2str(compStruct.subsystems.matrix(i,j));
        end
    end
    subsystem_comparison = str;
    subsystem_comparison
    fprintf('\n\n');
end
if plotResults==true
    plottingData = compStruct.subsystems.matrix./median(compStruct.subsystems.matrix,2);
    color_map = redblue(length(0:.01:2));
    h = genHeatMap(plottingData,compStruct.modelIDs,compStruct.subsystems.ID,'both','euclidean',color_map,[0,2]);
end

field = 'rxns';
all_rxns = catModelElements(models,field);
% Create binary matrix of reactions
binary_matrix = zeros(length(all_rxns),numel(models));
for i=1:numel(models)
    binary_matrix(:,i) = ismember(all_rxns,models{i}.rxns);
end
compStruct.structComp = squareform(1-pdist(binary_matrix','hamming'));
for i = 1:101
    color_map(i,:) = [(i-1)/100 0 0];
end
if plotResults == true
    h = genHeatMap(compStruct.structComp,compStruct.modelIDs,compStruct.modelIDs,'both','hamming',color_map,[0,1]);
end

end

%% Additional Functions
function A=getElements(models,field)
    A={};
    for i=1:numel(models)
       if isfield(models{i},field)
           A=[A;{models{i}.(field)}];
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

function B = compareSubsystems(models,field)
    % Compares number of reactions in each subsystem included in the models
    allSubsys = catModelElements(models,field);
    B = zeros(length(allSubsys),numel(models));
    for i=1:numel(models)
        for j=1:length(allSubsys)
            B(j,i) = sum(ismember(models{i}.subSystems,allSubsys(j)));
        end
    end
end

function C = catModelElements(models,field)
    % Creates a vector of all IDs included in models
    A = getElements(models,field);
    allElements = [];
    for i=1:numel(models)
        allElements = [allElements;A{i}];
    end
    C = unique(allElements);
end

function h = genHeatMap(data,colnames,rownames,clust_dim,clust_dist,col_map,col_bounds)
%GENHEATMAP  Generate a heatmap for a given matrix of data.
%
%   GENHEATMAP(DATA,COLNAMES,ROWNAMES,CLUSTER,COL_MAP,COL_BOUNDS) takes as
%   input a matrix of values and produces a heatmap using the built-in
%   PCOLOR function.
%
%--------------------------------- INPUTS ---------------------------------
%
% data        Numerical matrix.
%
% colnames    Names of DATA columns.
% 
% rownames    Names of DATA rows.
%
% clust_dim   'none' - the data will be plotted as provided (DEFAULT)
%             'rows' - cluster/rearrange the rows based on distance
%             'cols' - cluster/rearrange the columns based on distance
%             'both' - cluster/rearrange rows and columns based on distance
%
% clust_dist  Distance metric to be used for clustering, ignored if
%             CLUST_DIM is 'none'. Options are the same as those for
%             distance in, e.g., PDIST ('euclidean', 'hamming', etc.).
%             (DEFAULT = 'euclidean')
%
% col_map     Colormap, provided as string (e.g., 'parula', 'hot', 'etc.')
%             or an Nx3 RGB matrix of N colors.
%             (DEFAULT = 'parula')
%
% col_bounds  A 2-element vector with min and max values, to manually set
%             the bounds of the colormap.
%             (DEFAULT = min/max of DATA).
%
%
% Jonathan Robinson, 2018-03-06
%


% ***** additional adjustable parameters *****
grid_color = 'none';  % color of gridlines (RGB vector or string). e.g., [1 0 0] or 'r' (for red)
linkage_method = 'average';  % linkage algorithm for hierarchical clustering (see linkage function for more options)
% ********************************************


% handle input arguments
if nargin < 7 || isempty(col_bounds)
    col_bounds = [min(data(:)),max(data(:))];
    if nargin < 6 || isempty(col_map)
        col_map = 'parula';
        if nargin < 5 || isempty(clust_dist)
            clust_dist = 'euclidean';
            if nargin < 4 || isempty(clust_dim)
                clust_dim = 'none';
            elseif ~ismember(clust_dim,{'none','rows','cols','both'})
                error('%s is not a valid CLUST_DIM option. Choose "none", "rows", "cols", or "both".',clust_dim);
            end
        end
    end
end


% perform hierarchical clustering to sort rows (if specified)
if ismember(clust_dim,{'rows','both'})
    L = linkage(data,linkage_method,clust_dist);
    row_ind = optimalleaforder(L,pdist(data,clust_dist));
else
    row_ind = 1:size(data,1);
end
% perform hierarchical clustering to sort columns (if specified)
if ismember(clust_dim,{'cols','both'})
    L = linkage(data',linkage_method,clust_dist);
    col_ind = optimalleaforder(L,pdist(data',clust_dist));
else
    col_ind = 1:size(data,2);
end

% reorder data matrix according to clustering results
sortdata = data(row_ind,col_ind);
sortrows = rownames(row_ind);
sortcols = colnames(col_ind);

% check if data is square matrix with identical row and column names
if (length(colnames) == length(rownames)) && all(strcmp(colnames,rownames))
    % flip data so the diagonal is from upper left to lower right
    sortdata = fliplr(sortdata);
    sortcols = flipud(sortcols);
end

% pad data matrix with zeros (pcolor cuts off last row and column)
sortdata(end+1,end+1) = 0;

% generate pcolor plot
warning('off','MATLAB:hg:willberemoved');  % suppress warnings about tick label adjustment not being supported in future Matlab releases.

figure(); % Generate new figure?
a = axes;
set(a,'YAxisLocation','Right','XTick',[],'YTick', (1:size(sortdata,1))+0.5,'YTickLabels',sortrows);
set(a,'TickLength',[0 0],'XLim',[1 size(sortdata,2)],'YLim',[1 size(sortdata,1)]);
hold on

h = pcolor(sortdata);
set(h,'EdgeColor',grid_color);
set(gca,'XTick', (1:size(sortdata,2))+0.5);
set(gca,'YTick', (1:size(sortdata,1))+0.5);
set(gca,'XTickLabels',sortcols,'YTickLabels',sortrows);
set(gca,'XTickLabelRotation',90);
colormap(col_map);

if ~isempty(col_bounds)
    caxis(col_bounds);
end

xl = get(gca,'XLim');
yl = get(gca,'YLim');
plot(xl([1,1,2,2,1]),yl([1,2,2,1,1]),'k');

warning('on','MATLAB:hg:willberemoved');  % turn the warning back on
end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];
end