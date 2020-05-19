function [model,MNXrefMets,MNXrefRxns] = buildMNXref(type,allIDs,mnxPath,version,saveAsMNXref)
%buildMNXref  Construct a reference structure from MNX database files.
%
%   type            string, specifying the type of reference structure
%                   'mets'  include only metabolite-related fields.
%                   'rxns'  included only reaction-related fields
%                   'both'  include fields related to both reactions and
%                           metabolites. (opt, default 'both')
%   allIDs          string whether all IDs from different databases should be
%                   provided for each metabolite or reaction.
%                   'all'   all IDs from the ...._xref files are included
%                   'chebi' same as 'all', but filtered for secondary ChEBI IDs
%                   'pref'  only preferred IDs from ...._prop files (one per
%                           metabolite / reaction) are included. (opt, default
%                           'chebi')
%   mnxPath         string of path where MetaNetX reference files (chem_xref,
%                   chem_prop, reac_xref, reac_prop and chebiSecondary) are
%                   stored. (opt, default to RAVENdir/external/metanetx)
%
%   version         string, specifying the version of MNXref (opt, default '3.2',
%                   the latest version as of 2020-05-14)
%
%   saveAsMNXref    string whether to save structures MNXrefMets and MNXrefRxns as a
%                   single .mat file in the current directory.
%                   'yes'   save structures as MNXref.mat
%                   'no'    MNXref.mat file will not be created (opt, default 'no')         
%
% Usage: [model,MNXrefMets,MNXrefRxns] = buildMNXref(type,allIDs,mnxPath,saveAsMNXref)
%
% Eduard Kerkhoven, 2018-07-25
% Cheng Wei Quan (Eiden), 2020-05-14

if nargin<2
    allIDs='chebi';
end

if nargin<3
    [ST, I]=dbstack('-completenames');
    mnxPath=fileparts(fileparts(fileparts(ST(I).file)));
    mnxPath=fullfile(mnxPath,'external','metanetx');
end

if nargin<4
    version = '3.2';
end

if nargin<5
    saveAsMNXref = 'no';
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
        'Delimiter','\t','NumHeaderLines',385);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx = readtable(fullfile(mnxPath,'reac_prop.txt'),opts);
    
    opts=detectImportOptions(fullfile(mnxPath,'reac_xref.txt'),...
        'Delimiter','\t','NumHeaderLines',385);
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
        
    % extract source and ID information from MNX data
    if strcmp(allIDs,'pref')
        rxnSources = mnx.Source;
        rxnSources(~contains(rxnSources,':')) = {''};  % only keep entries with a colon, which indicates they have a source beyond MNX
        rxnSourceNames = regexprep(rxnSources,':.+$','');  % retrieve source name
        rxnSourceIDs = regexprep(rxnSources,'^.+:','');  % retrieve source ID
        mnxID2DBID = [mnx.MNX_ID,rxnSourceNames,rxnSourceIDs];
    end
    rxnSources = regexprep(mnx_xref.XREF,':.+$','');  % retrieve source name
    rxnSourceIDs = regexprep(mnx_xref.XREF,'^.+:','');  % retrieve source ID
    mnxID2extID = [mnx_xref.MNX_ID,rxnSources,rxnSourceIDs];
    
    % Remove any BIGG id that ends with small letter: localization, which
    % is not mapped in this function
    rmIds = strcmp(mnxID2extID(:,2),'bigg') & ~cellfun(@isempty,regexp(mnxID2extID(:,3),'.*[a-z]$'));
    rmIds = logical(rmIds + strcmp(mnxID2extID(:,2),'bigg') & startsWith(mnxID2extID(:,3),'R_'));
    % Only keep first entry for each reaction in rhea
    rheaIdx = find(strcmp(mnxID2extID(:,2),'rhea'));
    [~,keep,~] = unique(mnxID2extID(rheaIdx,1),'first');
    rheaRemove = ~ismember(rheaIdx,rheaIdx(keep));
    rmIds(rheaIdx(rheaRemove)) = 1;
    
    mnxID2extID(rmIds,:)=[];
    mnx_xref(rmIds,:)=[];
        
    sNames={'bigg','kegg','metacyc','reactome','rhea','sabiork','seed'};
    fNames={'rxnBIGGID','rxnKEGGID','rxnMetaCycID','rxnREACTOMEID',...
        'rxnRheaID','rxnSABIORKID','rxnSEEDID'};
    if strcmp(allIDs,'pref')
        mnxIDdb=mnxID2DBID;
    else
        mnxIDdb=mnxID2extID;
    end   
    for i=1:length(sNames)
        fprintf(' done.\nRetrieving %s reaction IDs...',upper(sNames{i}));
        currentDbOnly=ismember(mnxIDdb(:,2),sNames{i});
        currentmnxID2extID=mnxIDdb(currentDbOnly,:);
        
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
    model.rxnMetaNetXID = model.rxns;
    [~,ind] = ismember(mnxID2extID(:,2),sNames);
    mnxID2extID(:,2) = fNames(ind);  % rename to match field names
    model.rxnMetaNetXID2extID = mnxID2extID;
    fprintf(' done.\n');
    
end

%% Load metabolite properties
if ismember(type,{'mets','both'})
    fprintf('Loading metabolite data... ');
    opts=detectImportOptions(fullfile(mnxPath,'chem_prop.txt'),...
        'Delimiter','\t','NumHeaderLines',385);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx = readtable(fullfile(mnxPath,'chem_prop.txt'),opts);  % ~45 sec load time
    opts=detectImportOptions(fullfile(mnxPath,'chem_xref.txt'),...
        'Delimiter','\t','NumHeaderLines',385);
    opts.VariableNames(1)=regexprep(opts.VariableNames(1),'^x_','');
    mnx_xref = readtable(fullfile(mnxPath,'chem_xref.txt'),opts);
    fprintf('done.\n');

    fprintf('Processing metabolite data...');
    model.mets = mnx.MNX_ID;
    model.metMetaNetXID = mnx.MNX_ID;
    
    % extract information and add to model
    model.metNames = mnx.Description;
    model.metFormulas = mnx.Formula;
    model.metCharges = mnx.Charge;
    model.inchis = mnx.InChI;
    model.metSMILES = mnx.SMILES;
    
    % remove rows of mnx_xref that don't contain external IDs
    mnx_xref(~contains(mnx_xref.XREF,':') | contains(mnx_xref.XREF,'deprecated:'),:) = [];
    
    % extract source and ID information from MNX data
    if strcmp(allIDs,'pref')
        metSources = mnx.Source;
        metSources(~contains(metSources,':')) = {''};  % only keep entries with a colon, which indicates they have a source beyond MNX
        metSourceNames = regexprep(metSources,':.+$','');  % retrieve source name
        metSourceIDs = regexprep(metSources,'^.+:','');  % retrieve source ID
        mnxID2DBID = [mnx.MNX_ID,metSourceNames,metSourceIDs];
    end
    
    metSourceNamesX = regexprep(mnx_xref.XREF,':.+$','');  % retrieve source name
    metSourceIDsX = regexprep(mnx_xref.XREF,'^.+:','');  % retrieve source ID
    mnxID2extID = [mnx_xref.MNX_ID,metSourceNamesX,metSourceIDsX];
    
    rmIds = strcmp(mnxID2extID(:,2),'bigg') & startsWith(mnxID2extID(:,3),'M_');
    rmIds(find(strcmp(mnxID2extID(:,2),'kegg') & startsWith(mnxID2extID(:,3),'D'))) = 1;
    rmIds(find(strcmp(mnxID2extID(:,2),'kegg') & startsWith(mnxID2extID(:,3),'G'))) = 1;
    
    mnxID2extID(rmIds,:)=[];
    mnx_xref(rmIds,:)=[];

    % Split met name field by delimiters '|' and '; '. This takes a very
    % long time, and should only be done if necessary.
    % change semicolon+space delimiters to vertical bar "|"
    regexprep(mnx_xref.Description,'; ','|');
    
    % break processing into 10 chunks
    chunk_ind = 1:round(length(mnx_xref.Description)/10):length(mnx_xref.Description);
    chunk_ind(end) = length(mnx_xref.Description)+1;
    
    %split met names
    mnxDescr = cell(size(mnx_xref,1),1);
    fprintf(' done.\nSplitting synonymous MNX metabolite names: ');
    for i = 1:length(chunk_ind)-1
        fprintf('%u%% ',round(chunk_ind(i)/length(mnx_xref.Description)*100)+10);
        mnxDescr(chunk_ind(i):chunk_ind(i+1)-1) = cellfun(@(x) strsplit(x,'|'),mnx_xref.Description(chunk_ind(i):chunk_ind(i+1)-1),'UniformOutput',false);
    end
    fprintf('100%% done.\n');
    % flatten nested cells to obtain an Nx2 matrix of MNXID-name pairs
    fprintf('Indexing metabolite names to their MNX ID...');
    mnxIDs = cellfun(@(id,descr) repmat({id},1,numel(descr)),mnx_xref.MNX_ID,mnxDescr,'UniformOutput',false);
    model.mnxID2name = [[mnxIDs{:}]',[mnxDescr{:}]'];
    
    % remove repeated ID-name pairs
    [~,numericPair] = ismember(lower(model.mnxID2name),unique(lower(model.mnxID2name)));  % make numeric for faster processing
    [~,uniq_ind] = unique(numericPair,'rows');
    model.mnxID2name = model.mnxID2name(uniq_ind,:);
        
    % add model fields corresponding to available source names
    sNames = {'bigg','chebi','envipath','hmdb','kegg','lipidmaps','metacyc','reactome','sabiork','seed','slm'};
    fNames = {'metBIGGID','metChEBIID','metEnviPathID','metHMDBID','metKEGGID','metLIPIDMAPSID','metMetaCycID','metREACTOMEID','metSABIORKID','metSEEDID','metSLMID'};
    if strcmp(allIDs,'pref')
        mnxIDdb=mnxID2DBID;
    else
        mnxIDdb=mnxID2extID;
    end    
    for i=1:length(sNames)
        fprintf(' done.\nRetrieving %s metabolite IDs...',upper(sNames{i}));
        currentDbOnly=ismember(mnxIDdb(:,2),sNames{i});
        currentmnxID2extID=mnxIDdb(currentDbOnly,:);
        
        if strcmp(sNames(i),'bigg')
            rmIds=startsWith(currentmnxID2extID(:,3),'M_');
            currentmnxID2extID(rmIds,:)=[];
        end
    
        if i==2 && strcmp(allIDs,'chebi')
            chebi=fileread([mnxPath,'\chebi.dat']);
            chebi=strsplit(chebi,'\n');
            secChebi=~ismember(currentmnxID2extID(:,3),chebi);
            currentmnxID2extID(secChebi,:)=[];
        end
        
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

%% Process fields with empty cells/nested cells in model
if ismember(type,{'rxns','both'})
[MNXrefRxns.BiGGMNXid,MNXrefRxns.BiGGxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnBIGGID);
[MNXrefRxns.KEGGMNXid,MNXrefRxns.KEGGxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnKEGGID);
[MNXrefRxns.MetaCycMNXid,MNXrefRxns.MetaCycxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnMetaCycID);
[MNXrefRxns.ReactomeMNXid,MNXrefRxns.Reactomexref] = extractNestedCell(model.rxnMetaNetXID,model.rxnREACTOMEID);
[MNXrefRxns.RheaMNXid,MNXrefRxns.Rheaxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnRheaID);
[MNXrefRxns.SABIORKMNXid,MNXrefRxns.SABIORKxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnSABIORKID);
[MNXrefRxns.SEEDMNXid,MNXrefRxns.SEEDxref] = extractNestedCell(model.rxnMetaNetXID,model.rxnSEEDID);
MNXrefRxns.Version = version;
end

if ismember(type,{'mets','both'})
[MNXrefMets.BiGGMNXid,MNXrefMets.BiGGxref] = extractNestedCell(model.metMetaNetXID,model.metBIGGID);
[MNXrefMets.ChEBIMNXid,MNXrefMets.ChEBIxref] = extractNestedCell(model.metMetaNetXID,model.metChEBIID);
[MNXrefMets.envipathMNXid,MNXrefMets.envipathxref] = extractNestedCell(model.metMetaNetXID,model.metEnviPathID);
[MNXrefMets.HMDBMNXid,MNXrefMets.HMDBxref] = extractNestedCell(model.metMetaNetXID,model.metHMDBID);
[MNXrefMets.KEGGMNXid,MNXrefMets.KEGGxref] = extractNestedCell(model.metMetaNetXID,model.metKEGGID);
[MNXrefMets.LIPIDMAPSMNXid,MNXrefMets.LIPIDMAPSxref] = extractNestedCell(model.metMetaNetXID,model.metLIPIDMAPSID);
[MNXrefMets.MetaCycMNXid,MNXrefMets.MetaCycxref] = extractNestedCell(model.metMetaNetXID,model.metMetaCycID);
[MNXrefMets.ReactomeMNXid,MNXrefMets.Reactomexref] = extractNestedCell(model.metMetaNetXID,model.metREACTOMEID);
[MNXrefMets.SABIORKMNXid,MNXrefMets.SABIORKxref] = extractNestedCell(model.metMetaNetXID,model.metSABIORKID);
[MNXrefMets.SEEDMNXid,MNXrefMets.SEEDxref] = extractNestedCell(model.metMetaNetXID,model.metSEEDID);
[MNXrefMets.SLMMNXid,MNXrefMets.SLMxref] = extractNestedCell(model.metMetaNetXID,model.metSLMID);
MNXrefMets.Version = version;
end

if strcmp(saveAsMNXref,'yes')
    save('MNXref','MNXrefMets','MNXrefRxns');
end
end

