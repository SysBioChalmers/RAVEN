function model = buildMNXref(type,mnxPath)
%buildMNXref  Construct a reference structure from MNX database files.
%
%   type        string, specifying the type of reference structure
%               'mets'  include only metabolite-related fields.
%               'rxns'  included only reaction-related fields
%               'both'  include fields related to both reactions and
%                       metabolites. (opt, default 'both')
%   mnxPath     string of path where MetaNetX reference files (chem_xref,
%               chem_prop, reac_xref and reac_prop) are stored. (opt,
%               default to RAVENdir/external/metanetx)
%
% Usage: model = buildMNXref(model_type,mnxPath)
%
% Eduard Kerkhoven, 2018-07-22

if ~exist('mnxPath','var')
    [ST, I]=dbstack('-completenames');
    mnxPath=fileparts(fileparts(fileparts(ST(I).file)));
    mnxPath=fullfile(mnxPath,'external','metanetx');
end

files={'chem_xref','chem_prop','reac_prop','reac_xref'};
if strcmp(type,'mets')
    files=files(1:2);
elseif strcmp(type,'rxns')
    files=files(3:4);
end
for i=1:length(files)
    if exist(fullfile(mnxPath,[files{i},'.tsv']),'file')
        movefile(fullfile(mnxPath,[files{i},'.tsv']),...
            fullfile(mnxPath,[files{i},'.txt']));
    end
    if ~exist(fullfile(mnxPath,[files{i},'.txt']), 'file')
        fprintf('File %s cannot be found and will be downloaded from MetaNetX.org...\n',files{i});
        websave(fullfile(mnxPath,[files{i},'.txt']),...
            ['https://www.metanetx.org/cgi-bin/mnxget/mnxref/',files{i},'.tsv']);
    end
end

% handle inputs
if nargin < 1
    type = 'both';
end

%% Initialize model
model=[];

%% Load reaction properties
if ismember(type,{'rxns','both'})
    fprintf('Loading reaction data files...');
    opts=detectImportOptions(fullfile(mnxPath,'reac_prop.txt'),...
        'Delimiter','\t','NumHeaderLines',365);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx = readtable(fullfile(mnxPath,'reac_prop.txt'),opts);
    
    opts=detectImportOptions(fullfile(mnxPath,'reac_xref.txt'),...
        'Delimiter','\t','NumHeaderLines',365);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx_xref = readtable(fullfile(mnxPath,'reac_xref.txt'),opts);
    
    % add a field that contains a list of mets that participate in each reaction
    rxnMets = cellfun(@(r) regexp(r,'(\w+)@\w+','tokens'),mnx.Equation,'UniformOutput',false);
    model.rxnMets = cellfun(@(r) unique([r{:}]),rxnMets,'UniformOutput',false);
    
    % extract reaction MNX IDs and other information
    model.rxns = mnx.MNX_ID;
    model.rxnBalanced = mnx.Balance;
    model.eccodes = mnx.EC;
    model.rxnEqnMNX = mnx.Equation;
    model.rxnEqnNames = mnx.Description;
    
    % get other rxn IDs
    %remove deprecated and other lines
    mnx_xref(~contains(mnx_xref.XREF,':') | contains(mnx_xref.XREF,'deprecated:'),:) = [];
    
    % extract source and ID information from MNX xref data
    rxnSources = regexprep(mnx_xref.XREF,':.+$','');  % retrieve source name
    rxnSourceIDs = regexprep(mnx_xref.XREF,'^.+:','');  % retrieve source ID
    mnxID2extID = [mnx_xref.MNX_ID,rxnSources,rxnSourceIDs];
    
    sNames={'bigg','kegg','metacyc','reactome','rhea','sabiork','seed'};
    fNames={'rxnBIGGID','rxnKEGGID','rxnMetaCycID','rxnREACTOMEID',...
        'rxnRheaID','rxnSABIORKID','rxnSEEDID'};
    
    for i=1:length(sNames)
        fprintf(' done.\nRetrieving %s reaction IDs...',upper(sNames{i}));
        currentDbOnly=ismember(mnxID2extID(:,2),sNames{i});
        currentmnxID2extID=mnxID2extID(currentDbOnly,:);
        
        [mnxs,indices,indices2] = unique(currentmnxID2extID(:,1),'stable');
        counts = hist(indices2, 1:size(indices));
        mnxs(counts==1,2)=currentmnxID2extID(indices(counts==1),3);
        multiMnxs=find(counts>1);
        
        for j=1:length(multiMnxs)
            idx=strcmp(mnxs(multiMnxs(j)),currentmnxID2extID(:,1));
            mnxs(multiMnxs(j),2)={currentmnxID2extID(idx,3)};
        end
        
        [~,idx]=ismember(mnxs(:,1),model.rxns);
        DBstruct=cell(length(model.rxns),1);
        DBstruct(idx)=mnxs(:,2);
        model.(fNames{i})=DBstruct;
    end
    model.rxnMNXID = model.rxns;
    fprintf(' done.\n');
end

%% Load metabolite properties
if ismember(type,{'mets','both'})
    fprintf('Loading metabolite data... ');
    opts=detectImportOptions(fullfile(mnxPath,'chem_prop.txt'),...
        'Delimiter','\t','NumHeaderLines',365);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx = readtable(fullfile(mnxPath,'chem_prop.txt'),opts);  % ~45 sec load time
    opts=detectImportOptions(fullfile(mnxPath,'chem_xref.txt'),...
        'Delimiter','\t','NumHeaderLines',365);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx_xref = readtable(fullfile(mnxPath,'chem_xref.txt'),opts);
    fprintf('done.\n');

    fprintf('Processing metabolite data...');
    model.mets = mnx.MNX_ID;
    model.metMNXID = mnx.MNX_ID;
    
    % extract information and add to model
    model.metNames = mnx.Description;
    model.metFormulas = mnx.Formula;
    model.metCharges = mnx.Charge;
    model.inchis = mnx.InChI;
    model.metSMILES = mnx.SMILES;
    
    % get other IDs
    metSources = mnx.Source;
    metSources(~contains(metSources,':')) = {''};  % only keep entries with a colon, which indicates they have a source beyond MNX
    metSourceNames = regexprep(metSources,':.+$','');  % retrieve source name
    metSourceIDs = regexprep(metSources,'^.+:','');  % retrieve source ID
    
    % remove rows of mnx_xref that don't contain external IDs
    mnx_xref(~contains(mnx_xref.XREF,':') | contains(mnx_xref.XREF,'deprecated:'),:) = [];
    
    % extract source and ID information from MNX xref data
    metSourceNamesX = regexprep(mnx_xref.XREF,':.+$','');  % retrieve source name
    metSourceIDsX = regexprep(mnx_xref.XREF,'^.+:','');  % retrieve source ID
    mnxID2extID = [mnx_xref.MNX_ID,metSourceNamesX,metSourceIDsX];
    
    % Split met name field by delimiters '|' and '; '. This takes a very
    % long time, and should only be done if necessary.
    % change semicolon+space delimiters to vertical bar "|"
    regexprep(mnx_xref.Description,'; ','|');
    
    % break processing into 10 chunks
    chunk_ind = 1:round(length(mnx_xref.Description)/10):length(mnx_xref.Description);
    chunk_ind(end) = length(mnx_xref.Description)+1;
    
    % split met names
    mnxDescr = cell(size(mnx_xref,1),1);
    fprintf(' done.\nSplitting MNX metabolite names: ');
    for i = 1:length(chunk_ind)-1
        fprintf('%u%% ',round(chunk_ind(i)/length(mnx_xref.Description)*100)+10);
        mnxDescr(chunk_ind(i):chunk_ind(i+1)-1) = cellfun(@(x) strsplit(x,'|'),mnx_xref.Description(chunk_ind(i):chunk_ind(i+1)-1),'UniformOutput',false);
    end
    fprintf('100%% done\n');
    % flatten nested cells to obtain an Nx2 matrix of MNXID-name pairs
    fprintf('Organizing MNXID-name pairs... ');
    mnxIDs = cellfun(@(id,descr) repmat({id},1,numel(descr)),mnx_xref.MNX_ID,mnxDescr,'UniformOutput',false);
    model.mnxID2name = [[mnxIDs{:}]',[mnxDescr{:}]'];
    
    % remove repeated ID-name pairs
    [~,numericPair] = ismember(lower(model.mnxID2name),unique(lower(model.mnxID2name)));  % make numeric for faster processing
    [~,uniq_ind] = unique(numericPair,'rows');
    model.mnxID2name = model.mnxID2name(uniq_ind,:);
    
    % instead of splitting add unsplit name field as additional column
    %mnxID2extID = [mnxID2extID,mnx_xref.Description];
    
    % add model fields corresponding to available source names
    sNames = {'bigg','chebi','envipath','hmdb','kegg','lipidmaps','metacyc','reactome','sabiork','seed','slm'};
    fNames = {'metBiGGID','metChEBIID','metEnviPathID','metHMDBID','metKEGGID','metLIPIDMAPSID','metMetaCycID','metREACTOMEID','metSABIORKID','metSEEDID','metSLMID'};
    for i=1:length(sNames)
        fprintf(' done\nRetrieving %s reaction IDs...',upper(sNames{i}));
        currentDbOnly=ismember(mnxID2extID(:,2),sNames{i});
        currentmnxID2extID=mnxID2extID(currentDbOnly,:);
        
        [mnxs,indices,indices2] = unique(currentmnxID2extID(:,1),'stable');
        counts = hist(indices2, 1:size(indices));
        mnxs(counts==1,2)=currentmnxID2extID(indices(counts==1),3);
        multiMnxs=find(counts>1);
        
        for j=1:length(multiMnxs)
            idx=strcmp(mnxs(multiMnxs(j)),currentmnxID2extID(:,1));
            mnxs(multiMnxs(j),2)={currentmnxID2extID(idx,3)};
        end
        
        [~,idx]=ismember(mnxs(:,1),model.mets);
        DBstruct=cell(length(model.mets),1);
        DBstruct(idx)=mnxs(:,2);
        model.(fNames{i})=DBstruct;
    end
    fprintf(' done.\n')
    [~,ind] = ismember(mnxID2extID(:,2),sNames);
    mnxID2extID(:,2) = fNames(ind);  % rename to match field names
    model.mnxID2extID = mnxID2extID;
end

%% Final model adjustments
if ismember(type,{'mets','both'})
    % remove NAs from metFormulas
    model.metFormulas(ismember(model.metFormulas,'NA')) = {''};
end
