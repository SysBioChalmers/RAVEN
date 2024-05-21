function compStruct = compareMultipleModels(models,printResults,plotResults,groupVector,funcCompare,taskFile)
% compareMultipleModels
%   Compares two or more condition-specific models generated from the same
%   base model using high-dimensional comparisons in the reaction-space.
%
%   models              cell array of two or more models
%   printResults        true if the results should be printed on the screen
%                       (optional, default false)
%   plotResults         true if the results should be plotted
%                       (optional, default false)
%   groupVector         numeric vector or cell array for grouping similar 
%                       models, i.e. by tissue (optional, default, all models
%                       ungrouped)
%   funcCompare         logical, should a functional comparison be run
%                       (optional,default, false)
%   taskFile            string containing the name of the task file to use
%                       for the functional comparison (should be an .xls or 
%                       .xlsx file, required for functional comparison)
%
%   compStruct          structure that contains the comparison results
%       modelIDs        cell array of model ids
%       reactions       substructure containing reaction information
%           matrix          binary matrix composed of reactions (rows) in
%                           each model (column). This matrix is used as the
%                           input for the model comparisons.
%           IDs             list of the reactions contained in the reaction
%                           matrix.
%       subsystems      substructure containing subsystem information
%           matrix          matrix with comparison of number of rxns per
%                           subsystem
%           ID              vector consisting of names of all subsystems
%       structComp      matrix with pairwise comparisons of model structure
%                       based on (1-Hamming distance) between models
%       structCompMap   matrix with 3D tSNE (or MDS) mapping of model
%                       structures based on Hamming distances
%       funcComp        substructure containing function comparison results
%           matrix          matrix with PASS / FAIL (1 / 0) values for each
%                           task
%           tasks           vector containing names of all tasks
%
% Usage: compStruct=compareMultipleModels(models,printResults,...
%                       plotResults,groupVector,funcCompare,taskFile);

%% Stats toolbox required
if ~(exist('mdscale.m','file') && exist('pdist.m','file') && exist('squareform.m','file') && exist('tsne.m','file'))
    error('The MATLAB Statistics and Machine Learning Toolbox is required for this function')
end

%% Set up input defaults
if nargin < 2 || isempty(printResults)
    printResults=false;
end
if nargin < 3 || isempty(plotResults)
    plotResults=false;
end
if nargin < 4
    groupVector = [];
elseif ~isnumeric(groupVector)
    % convert strings to numeric groups
    [groupNames,~,groupVector] = unique(groupVector);
else
    % generate group names for vector of numbers
    groupNames = arrayfun(@num2str,unique(groupVector),'UniformOutput',false);
end
if nargin < 5 || isempty(funcCompare)
    funcCompare = false;
end
if nargin < 6
    taskFile = [];
else
    taskFile=char(taskFile);
end
if numel(models) <= 1
    EM = 'Cannot compare only one model. Use printModelStats if you want a summary of a model';
    dispEM(EM);
end
if isempty(taskFile) && funcCompare
    EM = 'Cannot perform the functional comparison without a task file. Specify taskFile or set funcCompare to FALSE.';
    dispEM(EM);
end

%% Set up model ID structure
compStruct.modelIDs = {};
fprintf('\n Getting model IDs \n')    
for i = 1:numel(models)
    if ~ischar(models{i}.id)  % to deal with non-character IDs (cells, strings, etc)
        compStruct.modelIDs{i,1} = models{i}.id{1};
    else
        compStruct.modelIDs{i,1} = models{i}.id;
    end
end
fprintf('*** Done \n\n')


%% Flatten models' subSystems field
% Convert from cell array of cells to cell array of strings
% NOTE: this function currently only recognizes one subSystem per reaction;
%       additional subSystems will be ignored!
for i = 1:numel(models)
    cells = cellfun(@iscell,models{i}.subSystems);
    models{i}.subSystems(cells) = cellfun(@(s) s{1}, models{i}.subSystems(cells), 'UniformOutput', false);
end


%% Compare models structure & function based on high-dimensional methods
% Compare number of reactions in each subsystem in each model using a heatmap
field = 'subSystems';
fprintf('\n Comparing subsystem utilization \n')
if any(~cellfun(@(m) isfield(m,field),models))
    fprintf('\nWARNING: At least one model does not contain the field "subSystems". \n')
    fprintf('         Skipping subsystem comparison. \n\n')
else
    [id,compMat] = compareModelField(models,field);
    compStruct.subsystems.ID = id;
    compStruct.subsystems.matrix = compMat;
    fprintf('*** Done \n\n')
    
    if printResults
        % This could use some cleaning up
        fprintf('*** Comparison of reaction subsystem populations:\n\n');
        
        nrow = min([15,numel(compStruct.subsystems.ID)]);
        ncol = min([10,numel(compStruct.modelIDs)]);
        summaryArray = [{field}, compStruct.modelIDs(1:ncol)'];
        summaryArray = [summaryArray; [compStruct.subsystems.ID(1:nrow), ...
            arrayfun(@num2str,compStruct.subsystems.matrix(1:nrow,1:ncol),'UniformOutput',false)]];
        
        charArray = [];
        for i = 1:size(summaryArray,2)
            charArray = [charArray, char(strcat(summaryArray(:,i),{'   '}))];
        end
        disp(charArray);
        if numel(compStruct.subsystems.ID) > 15
            fprintf('...\n');
        end
        fprintf('\n\n');
    end
    
    if plotResults==true
        % Plot all subsystems
        figure;
        plottingData = (compStruct.subsystems.matrix - mean(compStruct.subsystems.matrix,2))./mean(compStruct.subsystems.matrix,2);
        color_map = redblue(length(0:.01:2));
        h = genHeatMap(plottingData',compStruct.subsystems.ID,compStruct.modelIDs,'both','euclidean',color_map,[-1,1]);
        title('Subsystem Coverage - all subsystems','FontSize',18,'FontWeight','bold')
        
        % Plot only subsystems with deviation from mean
        keepSubs = (sum(plottingData~=0,2) ~= 0);
        if sum(keepSubs) > 1
            figure;
            h_small = genHeatMap(plottingData(keepSubs,:)',compStruct.subsystems.ID(keepSubs),...
                compStruct.modelIDs,'both','euclidean',color_map,[-1,1]);
            title('Subsystem Coverage','FontSize',18,'FontWeight','bold')
            
            % Plot enrichment in subsystems with deviation from mean
            figure;
            color_map_bw = [1 1 1;0 0 0];
            h_enriched = genHeatMap(plottingData(keepSubs,:)',compStruct.subsystems.ID(keepSubs),...
                compStruct.modelIDs,'both','euclidean',color_map_bw,[-10^-4,10^-4]);
            title('Subsystem Enrichment','FontSize',18,'FontWeight','bold')
        end
    end
    
end

% Compare overall reaction structure across all models using a heatmap
field = 'rxns';
fprintf('\n Comparing model reaction correlations \n')

% Create binary matrix of reactions
[id,binary_matrix] = compareModelField(models,field);
compStruct.reactions.IDs = id;
compStruct.reactions.matrix = binary_matrix;

% calculate hamming similarity
compStruct.structComp = 1-squareform(pdist(binary_matrix','hamming'));
fprintf('*** Done \n\n')
if plotResults == true
    color_map = [ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
    figure;
    h = genHeatMap(compStruct.structComp,compStruct.modelIDs,compStruct.modelIDs,'both','euclidean',color_map);
    title('Structural Similarity','FontSize',18,'FontWeight','bold')
end

% Compare overall reaction structure across all models using tSNE projection
rng(42) % For consistency
fprintf('\n Comparing model reaction structures \n')
if exist('tsne') > 0
    proj_coords = tsne(double(binary_matrix'),'Distance','hamming','NumDimensions',3); % 3D
    compStruct.structCompMap = proj_coords;
    axis_labels = {'tSNE 1';'tSNE 2';'tSNE 3'};
else % Seems odd to use mdscale if tsne is not found, as both are
     % distributed with stats toolbox. If tsne is not present, then also
     % mdscale is missing. Anyway, will leave this for now.
    [proj_coords,stress,disparities] = mdscale(pdist(double(binary_matrix'),'hamming'),3);
    compStruct.structCompMap = proj_coords;
    axis_labels = {'MDS 1';'MDS 2';'MDS 3'};
end
fprintf('*** Done \n\n')

% plot structure comparison results
if plotResults == true
    figure; hold on;
    if ~isempty(groupVector)
        color_data = groupVector;
        if length(groupNames) <= 7
            % "lines" colormap only has 7 unique colors
            color_palette = lines(length(groupNames));
        else
            color_palette = parula(length(groupNames));
        end
        colormap(color_palette);
    else
        color_data = 'k';
    end
    scatter3(proj_coords(:,1),proj_coords(:,2),proj_coords(:,3),35,color_data,'filled');
    view(135,25);  % to make it obvious that it is a 3D plot
    xlabel(axis_labels{1}); ylabel(axis_labels{2}); zlabel(axis_labels{3});
    set(gca,'FontSize',14,'LineWidth',1.25);
    title('Structural Comparison','FontSize',18,'FontWeight','bold')
    
    % add legend
    if ~isempty(groupVector)
        for i = 1:length(groupNames)
            h(i) = scatter3([],[],[],35,color_palette(i,:),'filled');
        end
        legend(h,groupNames);
    end
end

% Compare model functions by assessing their capacity to perform tasks
if funcCompare == true && ~isempty(taskFile)
    fprintf('\n Checking model performance on specified tasks. \n')
    taskStructure=parseTaskList(taskFile);
    for i = 1:numel(models)
        fprintf('\n Checking model # %.0f \n',i)
        taskReport{i} = checkTasks(models{i},[],false,false,false,taskStructure);
    end    
    
    % Save results
    taskMatrix = zeros(length(taskReport{1}.ok),numel(taskReport));
        for i = 1:numel(taskReport)
            taskMatrix(:,i) = taskReport{i}.ok;
        end
    compStruct.funcComp.matrix = taskMatrix;
    compStruct.funcComp.tasks = taskReport{1}.description;
    fprintf('*** Done \n\n')
   
    % Plot results
    if plotResults == true
        figure;
        color_map_bw = [1 1 1;0 0 0];
        h_enriched = genHeatMap(taskMatrix,compStruct.modelIDs,...
            taskReport{1}.description,'both','euclidean',color_map_bw,[0,1]);
        title('Functional Comparison - All Tasks','FontSize',18,'FontWeight','bold')
        
        figure;
        color_map_bw = [1 1 1;0 0 0];
        h_enriched = genHeatMap(taskMatrix(intersect(find(sum(taskMatrix,2)~=numel(models)),find(sum(taskMatrix,2)~=0)),:),...
            compStruct.modelIDs,taskReport{1}.description(intersect(find(sum(taskMatrix,2)~=numel(models)),find(sum(taskMatrix,2)~=0))),...
            'both','euclidean',color_map_bw,[0,1]);
        title('Functional Similarity','FontSize',18,'FontWeight','bold')
    end
end

end

%% Additional Functions

function [id,compMat] = compareModelField(models,field)
    % Generates a list of unique field entries and a matrix quantifying the
    % number of appearances of each field entry in each model
    
    % get unique list of field entries
    hasfield = cellfun(@(m) isfield(m,field),models);
    id = cellfun(@(m) m.(field),models(hasfield),'UniformOutput',false);
    id = unique(vertcat(id{:}));
    
    % assemble matrix comparing frequency of each entry in each model
    compMat = zeros(numel(id),numel(models));
    for i = 1:numel(models)
        [~,entryIndex] = ismember(models{i}.(field),id);  % get index of each field entry in the unique id list
        compMat(:,i) = histcounts(entryIndex, 0.5:1:(numel(id)+0.5));  % determine the frequency at which each index appears
    end
end


function h = genHeatMap(data,colnames,rownames,clust_dim,clust_dist,col_map,col_bounds,grid_color)
%genHeatMap  Generate a heatmap for a given matrix of data.
%
% Usage:
%
%   genHeatMap(data,colnames,rownames,clust_dim,clust_dist,col_map,col_bounds,grid_color);
%
% Inputs:
%
% data        Numerical matrix.
%
% colnames    Cell array of data column names.
% 
% rownames    Cell array of data row names.
%
% clust_dim   'none' - the data will be plotted as provided (DEFAULT)
%             'rows' - cluster/rearrange the rows based on distance
%             'cols' - cluster/rearrange the columns based on distance
%             'both' - cluster/rearrange rows and columns based on distance
%
% clust_dist  Distance metric to be used for clustering, ignored if
%             clust_dim is 'none'. Options are the same as those for
%             distance in, e.g., PDIST ('euclidean', 'hamming', etc.).
%             (DEFAULT = 'euclidean')
%
% col_map     Colormap, provided as string (e.g., 'parula', 'hot', 'etc.')
%             or an Nx3 RGB matrix of N colors.
%             (DEFAULT = 'hot')
%
% col_bounds  A 2-element vector with min and max values, to manually set
%             the bounds of the colormap.
%             (DEFAULT = min/max of data).
%
% grid_color  Color of the grid surrounding the heatmap cells.
%             (DEFAULT = 'none')
%
%

% handle input arguments
if nargin < 4 || isempty(clust_dim)
    clust_dim = 'none';
elseif ~ismember(clust_dim,{'none','rows','cols','both'})
    error('%s is not a valid CLUST_DIM option. Choose "none", "rows", "cols", or "both".',clust_dim);
end
if nargin < 5 || isempty(clust_dist)
    clust_dist = 'euclidean';
end
if nargin < 6 || isempty(col_map)
    col_map = 'hot';
end
if nargin < 7 || isempty(col_bounds)
    col_bounds = [min(data(:)),max(data(:))];
end
if nargin < 8
    grid_color = 'none';
end

% perform hierarchical clustering to sort rows (if specified)
linkage_method = 'average';
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