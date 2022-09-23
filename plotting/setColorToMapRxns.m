function [modifiedMap, xmlMap, model2] = setColorToMapRxns (mapFileName, model, v1, rxnsFluxTask, rxnLineWidth, rxnLineColor, molFillColor)
% USAGE:
% [modifiedMap, xmlMap, model2] =
%   setColorToMapRxns (mapFileName, model, v1, rxnsFluxTask, rxnLineWidth, rxnLineColor, molFillColor)
%
% INPUTS:
% mapFileName           relative path to XML map (e.g. 'yeastMap1.0.xml')
% model                 model structure
% v1                    vector of fluxes
% rxnsFluxTask          1 = Flux
%                       2 = Flux essentiality; map reactions causing any
%                           effect on growth (i.e. obj function)
%                       3 = Flux essentiality; map reactions causing no
%                           growth(i.e. obj function) when deleted, i.e.
%                           essential reactions
%                       4 = Flux essentiality; map reactions with not any
%                           effect on growth (i.e. obj function) when
%                           deleted
%                       5 = Omics
% ADDITIONAL INPUTS:
% rxnLineWidth          reaction line width; default = 5
% rxnLineColor          reaction line color; default = Pathways; write
%                       'Subsystems' to color five subsystems only (central,
%                       AA, lipids, cofactors, nuclotides and energy).
%                       'Subsystems' option is the deault option at task 5
% molFillColor          molecule fill color; default = true (i.e. white)
%                       molecule fill color; true (i.e. like CellDesigner default)
% OUTPUTS:
% modifiedMap           map with MATLAB structure, and modified reactions' line color and width
% xmlMap                xml structure of the original map obtained from mapFileName
% model2                model with modified subsystems

if ~(exist('transformXML2Map.m','file')==2)
    error('COBRA Toolbox is required, go to https://opencobra.github.io/cobratoolbox/')
end

if nargin < 5 || rxnLineWidth == 0
    rxnLineWidth = 0;
    rxnLineColor = '';
    molFillColor = true;
elseif nargin < 6
    rxnLineColor = '';
    molFillColor = true;
elseif nargin < 7
    molFillColor = true;
end
if rxnsFluxTask == 5
    rxnLineColor = 'Subsystems';
end

[xmlMap, mapMap] = transformXML2Map(mapFileName);
if molFillColor == true
    modifiedMap = unifyMetabolicMapCD(mapMap); % molColor whrite and rxnColor = ligthgray
elseif molFillColor == false
    modifiedMap = mapMap; % like CellDesigner default
end

%%%%%%%% Gettin reactions line width with normalized fluxes %%%%%%%%%%%%%
model2 = model; % to not overwrite model
flxs(:, 1) = v1; % to not overwrite v1
absFlxs = abs(flxs);
rxnLineWidthInMap = absFlxs / max(absFlxs); %normalize fluxes line width
if rxnsFluxTask == 1 && rxnLineWidth == 0
    rxnLineWidthInMap(rxnLineWidthInMap > 0.1 & rxnLineWidthInMap <=  1) = 50;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.01 & rxnLineWidthInMap <= 0.1) = 25;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.001 & rxnLineWidthInMap <= 0.01) = 16;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.0001 & rxnLineWidthInMap <= 0.001) = 10;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.00001 & rxnLineWidthInMap <= 0.0001) = 6;
elseif rxnsFluxTask == 2 && rxnLineWidth == 0
    rxnLineWidthInMap(rxnLineWidthInMap > 0.90 & rxnLineWidthInMap <=  0.95) =  6;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.80 & rxnLineWidthInMap <= 0.90) =  8;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.50 & rxnLineWidthInMap <= 0.80) = 11;
    rxnLineWidthInMap(rxnLineWidthInMap > 0.25 & rxnLineWidthInMap <= 0.50) = 15;
    rxnLineWidthInMap(rxnLineWidthInMap >= 0.00 & rxnLineWidthInMap <= 0.25) = 20;
elseif rxnsFluxTask == 1 && rxnLineWidth ~= 0.0
    rxnLineWidthInMap(rxnLineWidthInMap == 1.0) = rxnLineWidth;
elseif rxnsFluxTask == 2 && rxnLineWidth ~= 0.0
    rxnLineWidthInMap(rxnLineWidthInMap < 0.95) = rxnLineWidth;
elseif rxnsFluxTask == 3 && rxnLineWidth == 0.0
    rxnLineWidthInMap(rxnLineWidthInMap == 0.0) = 5;
elseif rxnsFluxTask == 3 && rxnLineWidth ~= 0.0
    rxnLineWidthInMap(rxnLineWidthInMap == 0.0) = rxnLineWidth;
elseif rxnsFluxTask == 4 && rxnLineWidth == 0.0
    rxnLineWidthInMap(rxnLineWidthInMap > 0.95 & rxnLineWidthInMap <= 1.0) = 5.0;
elseif rxnsFluxTask == 4 && rxnLineWidth ~= 0.0
    rxnLineWidthInMap(rxnLineWidthInMap > 0.95 & rxnLineWidthInMap <= 1.0) = rxnLineWidth;
elseif rxnsFluxTask == 5 && rxnLineWidth == 0.0
    cutValue = 2.0; % transcriptomics minimum cut value of log change
    rxnLineWidthInMap = absFlxs;
    rxnLineWidthInMap(rxnLineWidthInMap > 5.0) = 30;
    rxnLineWidthInMap(rxnLineWidthInMap > 4.0 & rxnLineWidthInMap <= 5.0) = 25;
    rxnLineWidthInMap(rxnLineWidthInMap > 3.0 & rxnLineWidthInMap <= 4.0) = 20;
    rxnLineWidthInMap(rxnLineWidthInMap > 2.0 & rxnLineWidthInMap <= 3.0) = 15;
    rxnLineWidthInMap(rxnLineWidthInMap > 1.0 & rxnLineWidthInMap <= 2.0) = 10;
    rxnLineWidthInMap(rxnLineWidthInMap >= 0.5 & rxnLineWidthInMap <= 1.0) =  5;
end
for i = 1:length(modifiedMap.rxnName)
    a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
    if isempty(a)
        modifiedMap.rxnWidth{i} = 1;
    else
        modifiedMap.rxnWidth{i} = rxnLineWidthInMap(a);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Cofactors not included in the analyses %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Compartments included in the analyses %%%%%%%%%%%%%%%%%%%%%%%%%

% propossed cofactors list for PART A
cofactorsListA ={'h'; 'h2o'; 'nad'; 'nadh'; 'nadp'; 'nadph';...
    'fad'; 'fadh'; 'ppi'; 'adp'; 'gtp'; 'gdp'; 'gmp'; ...
    'utp'; 'udp'; 'ump'; 'ctp'; 'cdp'; 'cmp'};

% propossed cofactors list for PART B
cofactorsListB ={'h'; 'h2o'; 'co2'; 'coa'; 'accoa'; 'nad'; 'nadh'; 'nadp'; 'nadph';...
    'fad'; 'fadh'; 'pi'; 'ppi'; 'atp'; 'adp'; 'amp'; 'gtp'; 'gdp'; 'gmp';...
    'utp'; 'udp'; 'ump'; 'ctp'; 'cdp'; 'cmp'; 'gthox'; 'gthrd_c'; 'amet_c';...
    'ahcys_c'; 'glu__L'; 'gln__L'; 'akg'; 'nh4'};

% compartments list
compartmentsList = {'_e'; '_c'; '_m'; '_n'; '_x'; '_v'; '_r'; '_g'};

% include compartment to cofactor lists A and B
cofactorsListWcompA = cell(length(cofactorsListA), 1);
c = 1; % general counter
for i = 1:length(cofactorsListA)
    for j = 1:length(compartmentsList)
        cofactorsListWcompA{c, 1} = horzcat(cofactorsListA{i},compartmentsList{j});
        c = c + 1;
    end
end
cofactorsListWcompB = cell(length(cofactorsListB), 1);
c = 1;
for i = 1:length(cofactorsListB)
    for j = 1:length(compartmentsList)
        cofactorsListWcompB{c, 1}= horzcat(cofactorsListB{i},compartmentsList{j});
        c = c + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% getting transport reactions
[nrxns,~] = size(model2.rxns);
transRxns = cell(nrxns, 1);
for j = 1:size(model2.rxns)
    if true(contains(model2.subSystems {j}, 'S_Transport'))
        transRxns{j, 1} = model2.rxns{j};
    end
end
transRxns = transRxns(~cellfun('isempty',transRxns));

% getting subStystems without transport
SubSystems = (unique(getModelSubSystems(model2)));
[nSubSystems, ~] = size(SubSystems);
subSystemsNoTrans = cell(nSubSystems, 1);
for i = 1:nSubSystems
    if ~contains(SubSystems{i,1},'S_Transport')&& ...
            ~contains(SubSystems{i, 1},'Extracellular exchange') &&...
            ~contains(SubSystems{i, 1},'S_tRNA_charging')
        subSystemsNoTrans{i, 1} = SubSystems{i, 1};
    else
        subSystemsNoTrans{i, 1} = {};
    end
end
subSystemsNoTrans = subSystemsNoTrans(~cellfun('isempty',subSystemsNoTrans));

% defining subSystems of trasnport reactions to be colored, these
% transport reactions complete reaction conections in subSystems which
% require conexions among diferents compartments
[nTransRxns,~] = size(transRxns);
[nSubSystemsNoTrans, ~] = size (subSystemsNoTrans);
totalRxnsToIncludeA = cell(nTransRxns,1);
for i = 1:nSubSystemsNoTrans
    subSystemRxns = findRxnsFromSubSystem(model2, subSystemsNoTrans(i));
    metsSubSystemRxns = findMetsFromRxns(model2, subSystemRxns);
    [nMetsSubSystemRxns, ~] = size (metsSubSystemRxns);
    % getting unique compartments
    metsCompartment = cell(nMetsSubSystemRxns,1);
    for j = 1:nMetsSubSystemRxns
        str = metsSubSystemRxns{j};
        [~, s] = size(str);
        metsCompartment{j,1} = str(s-1:s);
    end
    metsCompartment = (unique(metsCompartment));
    [nUniqueMetsComp, ~] = size(metsCompartment);
    
    if nUniqueMetsComp > 1
        % getting duplicated mets like 'pyr_c' and 'pyr_m' which are 'pyr'
        % and pyr after erasing the compartment signature
        metsFromSubSystemRxnsNoComp = cell(nMetsSubSystemRxns,1);
        for j = 1:nMetsSubSystemRxns
            str = metsSubSystemRxns{j};
            [~, s] = size(str);
            metsFromSubSystemRxnsNoComp{j, 1} = str(1:s-2); % erasing comp signature
        end
        [~, o, ~] = unique(metsFromSubSystemRxnsNoComp,'first');
        indexToDupes = find(not(ismember(1:numel(metsFromSubSystemRxnsNoComp),o)));
        [~, nIndexToDupes] = size(indexToDupes);
        duplicatedMets = cell(nIndexToDupes,1);
        for j = 1:nIndexToDupes
            duplicatedMets(j, 1) = metsFromSubSystemRxnsNoComp (indexToDupes(1, j));
        end
        
        % getting reactions from duplicated mets
        dulicatedRxnMets = cell(nIndexToDupes * nUniqueMetsComp,1);
        c = 1; % duplicatedRxnMets counter
        for j = 1:nIndexToDupes
            str = duplicatedMets{j, 1};
            if ~ismember(str, cofactorsListA)
                for k = 1:nUniqueMetsComp
                    str1 = metsCompartment{k};
                    if ~ismember(str1, compartmentsList(2)) % compartment '_c'
                        dulicatedRxnMets{c, 1} = {horzcat(str, '_c'), horzcat(str, str1)};
                        c = c + 1;
                    end
                end
            end
        end
        dulicatedRxnMets = dulicatedRxnMets(~cellfun('isempty',dulicatedRxnMets));
        [ndulicatedRxnMets, ~] = size(dulicatedRxnMets);
        % from duplicated reactions, get reactions to subscribe a reactions
        % subsystem from non transport subsystems
        rxnsToInclude = cell(ndulicatedRxnMets, 1);
        for j = 1:ndulicatedRxnMets
            metsToFindInTrasnportRxns = dulicatedRxnMets{j, 1};
            for k = 1:nTransRxns
                rxnFormula = printRxnFormula(model2, transRxns{k, 1}, false);
                if true(contains(string(rxnFormula), metsToFindInTrasnportRxns (1, 1)))
                    if true(contains(rxnFormula, metsToFindInTrasnportRxns (1, 2)))
                        rxnsToInclude{j, 1} = transRxns{k, 1};
                    end
                end
            end
        end
        rxnsToInclude = rxnsToInclude(~cellfun('isempty',rxnsToInclude));
        rxnsToInclude = unique(rxnsToInclude);
        [s, ~] = size(rxnsToInclude);
        totalRxnsToIncludeA = [totalRxnsToIncludeA; rxnsToInclude];
        
        for j = 1:s
            rxnID = findRxnIDs(model2, rxnsToInclude{j});
            model2.subSystems(rxnID) = {subSystemsNoTrans(i)};
        end
        
    end
end
totalRxnsToIncludeA = totalRxnsToIncludeA(~cellfun('isempty',totalRxnsToIncludeA));
totalRxnsToIncludeA = unique(totalRxnsToIncludeA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get subsystems category. This is based on the number of conexions that a
% subsystem includes in the genome scale metabolic model2
subsLength ={};
for i = 1:length(subSystemsNoTrans)
    rxnsSubs = findRxnsFromSubSystem(model2,subSystemsNoTrans(i));
    metsSubsRxns = findMetsFromRxns(model2, rxnsSubs);
    metsSubsRxnsNoCof = {};
    for j = 1:length(metsSubsRxns)
        if ~ismember(metsSubsRxns(j),cofactorsListWcompB)
            metsSubsRxnsNoCof{j,1} = char(metsSubsRxns(j));
        end
    end
    metsSubsRxnsNoCof = metsSubsRxnsNoCof(~cellfun('isempty',metsSubsRxnsNoCof));
    
    metsRxns = findRxnsFromMets(model2, metsSubsRxnsNoCof);
    metsRxnsSubs = {};
    for j = 1:length(metsRxns)
        rxnID = findRxnIDs(model2,metsRxns(j));
        metsRxnsSubs{j,1} = string(model2.subSystems(rxnID));
    end
    metsRxnsSubs = unique(string(metsRxnsSubs));
    subsLength{i,1} = length(metsRxnsSubs);
    subsLength{i,2} = subSystemsNoTrans{i};
end

% Get  remaining transport reacctions, namely, reactions not included in
% PART A

[nTransRxns, ~] = size(transRxns);
remainingTransRxns = cell(nTransRxns, 1);
for i = 1:nTransRxns
    if ~ismember(transRxns(i),totalRxnsToIncludeA)
        remainingTransRxns{i,1} = transRxns{i, 1};
    end
end
remainingTransRxns = remainingTransRxns(~cellfun('isempty',remainingTransRxns));

% Assign remaining transport reactions
nonIncludSubs = {'S_tRNA_charging', 'S_Other',...
    'Biomass and maintenance functions', };
totalRxnsToIncludeAB = totalRxnsToIncludeA;
remTransRxnsWnewSubs = {};
remTransRxnsSubs = {};
rxnsToInclude = cell(length(remainingTransRxns), 1);
subSToInclude = cell(length(remainingTransRxns), 1);
for i = 1:length(remainingTransRxns)
    metsRemTransRxns = findMetsFromRxns(model2, remainingTransRxns(i));
    metsRxns = {};
    for j = 1:length(metsRemTransRxns)
        if true(contains(metsRemTransRxns{j},'_c') &&...
                ~ismember(metsRemTransRxns(j), cofactorsListWcompB))
            metsRxns = findRxnsFromMets(model2,metsRemTransRxns(j));
        end
    end
    metsRxns = metsRxns(~cellfun('isempty',metsRxns));
    
    rxnSubs = {};
    for j = 1:length(metsRxns)
        rxnID = findRxnIDs(model2, metsRxns(j));
        rxnSubs{j,1} = model2.subSystems{rxnID};
    end
    rxnSubs = rxnSubs(~cellfun('isempty',rxnSubs));
    
    r = [];
    rankMax = 0;
    subSToUseStr = cell(1,1);
    for j = 1:length(rxnSubs)
        if true(~ismember(string(metsRxns(j)),totalRxnsToIncludeA))
            for k = 1:length(subsLength)
                if true(and(ismember(rxnSubs{j}, subsLength{k,2}),...
                        ~ismember(rxnSubs{j}, nonIncludSubs)))
                    r(k, 1) = subsLength{k, 1};
                    subSToUseStr{1, 1} = string(subsLength(k,2));
                    if k > 1
                        if true(rankMax < r(k, 1))
                            rankMax = r(k, 1);
                            subSToUseStr{1, 1} = string(subsLength(k,2));
                        end
                    end
                end
            end
        end
    end
    if ~isempty(subSToUseStr(1,1))
        rxnsToInclude{i} = remainingTransRxns{i};
        subSToInclude{i} = char(subSToUseStr{1, 1});
        remTransRxnsWnewSubs{i,1} = remainingTransRxns(i);
        remTransRxnsSubs{i,1} = subSToUseStr(1, 1);
        totalRxnsToIncludeAB = [totalRxnsToIncludeAB; rxnsToInclude];
    end
end
totalRxnsToIncludeAB = totalRxnsToIncludeAB(~cellfun('isempty',totalRxnsToIncludeAB));
totalRxnsToIncludeAB = unique(totalRxnsToIncludeAB);

for j = 1:length(rxnsToInclude)
    rxnID = findRxnIDs(model2, rxnsToInclude{j});
    model2.subSystems(rxnID) = {subSToInclude(j)};
end

%%%%%%%%%%%%%%%%%%%%%%%% ASIGN SUBSYSTEMS' COLOR %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of reactions to include in the map
if rxnsFluxTask == 1
    for i = 1:length(model2.rxns)
        if absFlxs(i)==0 % zero flux reactions are not included
            model2.subSystems(i) = {{'NonIncluded'}};  % included at the end
        end
    end
elseif rxnsFluxTask == 2
    for i = 1:length(model2.rxns)
        if absFlxs(i) > 0.95 * max(absFlxs)% reactions causing growth (i.e. obj function) defect when deleting
            model2.subSystems(i) = {{'NonIncluded'}};   % included at the end
        end
    end
elseif rxnsFluxTask == 3
    for i = 1:length(model2.rxns)
        if absFlxs(i) > 0.0 % zero flux reactions only, i.e. essential reations
            model2.subSystems(i) = {{'NonIncluded'}};   % included at the end
        end
    end
elseif rxnsFluxTask == 4
    for i = 1:length(model2.rxns)
        if absFlxs(i) < 0.95 * max(absFlxs)% reactions fluxes causing no growth (i.e. obj function) defect when deleting
            model2.subSystems(i) = {{'NonIncluded'}};   % included at the end
        end
    end
elseif rxnsFluxTask == 5
    for i = 1:length(model2.rxns)
        if absFlxs(i) < cutValue % zero flux reactions are not included
            model2.subSystems(i) = {{'NonIncluded'}};  % included at the end
        end
    end
end

%See colors letter-code at:
%la.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/24497/versions/2/screenshot.PNG
%[xmlMap, mapMap] = transformXML2Map('YeastMapUnified1.0b.xml');
centralPathsList = {'S_GlycolysisGluconeogenesis', 'S_Pentose_Phosphate_Pathway',...
    'S_Pyruvate_Metabolism', 'S_Citric_Acid_Cycle', 'S_Taurine_Metabolism',...
    'S_Starch_and_Sucrose_Metabolism', 'S_Xylose_Metabolism',...
    'S_Fructose_and_Mannose_Metabolism', 'S_Galactose_metabolism',...
    'S_Arabinose_Metabolism', 'S_Alternate_Carbon_Metabolism', 'S_Methane_Metabolism'};
aminoacPathsList = {'S_Glutamate_metabolism', 'S_Glutamine_Metabolism',...
    'S_Alanine_and_Aspartate_Metabolism', 'S_Histidine_Metabolism',...
    'S_Threonine_and_Lysine_Metabolism', 'S_Glycine_and_Serine_Metabolism',...
    'S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism',...
    'S_Arginine_and_Proline_Metabolism', 'S_Cysteine_Metabolism',...
    'S_Methionine_Metabolism', 'S_Valine__Leucine__and_Isoleucine_Metabolism',...
    'S_Other_Amino_Acid_Metabolism', 'S_Asparagine_metabolism', 'S_Nitrogen_Metabolism'};
cofactorsPathsList = {'S_NAD_Biosynthesis', 'S_Pantothenate_and_CoA_Biosynthesis',...
    'S_Folate_Metabolism', 'S_Riboflavin_Metabolism', 'S_Quinone_Biosynthesis',...
    'S_Porphyrin_and_Chlorophyll_Metabolism'};
lipidsPathsList = {'S_Fatty_Acid__Biosynthesis', 'S_Fatty_Acid_Degradation',...
    'S_Fatty_Acid_Metabolism', 'S_Glycerolipid_Metabolism', 'S_Glycoprotein_Metabolism',...
    'S_Phospholipid_Biosynthesis', 'S_Phospholipid_Metabolism',...
    'S_Sphingolipid_Metabolism', 'S_Sterol_Metabolism'};
nucleotidsPathsList = {'S_Purine_and_Pyrimidine_Biosynthesis',...
    'S_Nucleotide_Salvage_Pathway', 'S_Thiamine_Metabolism', 'S_Pyridoxine_Metabolism'};
enrgyPathsList = {'Biomass and maintenance functions', 'S_Anaplerotic_reactions',...
    'S_Oxidative_Phosphorylation'};
colors = createColorsMap;

if isempty(rxnLineColor)||rxnLineColor == "Pathways"
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_GlycolysisGluconeogenesis', 'DARKBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Pentose_Phosphate_Pathway', 'STEELBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Pyruvate_Metabolism', 'DODGERBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Citric_Acid_Cycle', 'MEDIUMBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Taurine_Metabolism', 'MEDIUMBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Starch_and_Sucrose_Metabolism', 'BLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Xylose_Metabolism', 'STEELBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Fructose_and_Mannose_Metabolism', 'DODGERBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Galactose_metabolism', 'MEDIUMBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Arabinose_Metabolism', 'BLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Alternate_Carbon_Metabolism', 'DARKBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Methane_Metabolism', 'DARKBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Taurine_Metabolism', 'DARKBLUE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Anaplerotic_reactions', 'CRIMSON');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Complex_Alcohol_Metabolism', 'MAROON');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Oxidative_Phosphorylation', 'RED');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Purine_and_Pyrimidine_Biosynthesis', 'PINK');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Nucleotide_Salvage_Pathway', 'HOTPINK');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Thiamine_Metabolism', 'MAGENTA');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Pyridoxine_Metabolism', 'PINK');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_NAD_Biosynthesis', 'PLUM');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Pantothenate_and_CoA_Biosynthesis', 'VIOLET');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Folate_Metabolism', 'DARKORCHID');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Riboflavin_Metabolism', 'PURPLE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Quinone_Biosynthesis', 'MEDIUMORCHID');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Porphyrin_and_Chlorophyll_Metabolism', 'PLUM');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Glutamate_metabolism', 'ORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Glutamine_Metabolism', 'DARKORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Alanine_and_Aspartate_Metabolism', 'CORAL');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Histidine_Metabolism', 'GOLD');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Threonine_and_Lysine_Metabolism', 'ORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Glycine_and_Serine_Metabolism', 'DARKORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Tyrosine__Tryptophan__and_Phenylalanine_Metabolism', 'CORAL');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Arginine_and_Proline_Metabolism', 'GOLD');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Cysteine_Metabolism', 'ORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Methionine_Metabolism', 'DARKORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Valine__Leucine__and_Isoleucine_Metabolism', 'CORAL');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Other_Amino_Acid_Metabolism', 'GOLD');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Asparagine_metabolism', 'ORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Nitrogen_Metabolism', 'DARKORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Fatty_Acid__Biosynthesis', 'FORESTGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Fatty_Acid_Degradation', 'MEDIUMSEAGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Fatty_Acid_Metabolism', 'DARKGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Glycerolipid_Metabolism', 'SEAGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Glycoprotein_Metabolism', 'FORESTGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Phospholipid_Biosynthesis', 'MEDIUMSEAGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Phospholipid_Metabolism', 'DARKGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Sphingolipid_Metabolism', 'SEAGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Sterol_Metabolism', 'FORESTGREEN');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Other', 'BLACK');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'Biomass and maintenance functions', 'RED');
elseif rxnLineColor == "Subsystems"
    modifiedMap = colorSubsystem(modifiedMap, model2, centralPathsList, 'BLUE'); %mapMap to fluxMap
    modifiedMap = colorSubsystem(modifiedMap, model2, aminoacPathsList, 'ORANGE');
    modifiedMap = colorSubsystem(modifiedMap, model2, cofactorsPathsList, 'ORCHID');
    modifiedMap = colorSubsystem(modifiedMap, model2, lipidsPathsList, 'LIME');
    modifiedMap = colorSubsystem(modifiedMap, model2, nucleotidsPathsList, 'MAGENTA');
    modifiedMap = colorSubsystem(modifiedMap, model2, enrgyPathsList, 'RED');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Complex_Alcohol_Metabolism', 'GOLDENROD');
    modifiedMap = colorSubsystem(modifiedMap, model2, 'S_Other', 'GOLDENROD');
elseif ~(rxnLineColor == "Subsystems" || rxnLineColor == "Payhways" || isempty(rxnLineColor))
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        if ~isempty(a) && string(model2.subSystems(a)) ~= "NonIncluded"
            modifiedMap.rxnColor(i, 1) = {colors(rxnLineColor)};
        end
    end
end
if rxnLineColor == "Subsystems" || rxnLineColor == "Payhways" || isempty(rxnLineColor)
    h2o = {'H2Ot', 'H2Otm',	'H2Otn', 'H2Otv', 'H2Otp', 'H2Oter'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(h2o)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(h2o(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('CYAN')};
            end
        end
    end
    o2 = {'O2t', 'O2tm', 'O2ter', 'O2tp'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(o2)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(o2(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('CYAN')};
            end
        end
    end
    Pi = {'PIt2r', 'PIt2m', 'PIt5m', 'PIt2n', 'PIt2v', 'PIt2p', 'ATPtm_H', 'ATPtp_H', ...
        'ATP2tp_H', 'ATPS3g', 'ATPS3v', 'UTPtm', 'CTPtm', 'GTPt2m', 'ATPS'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(Pi)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(Pi(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('RED')};
            end
        end
    end
    co2 = {'CO2t', 'CO2tg', 'CO2tm', 'CO2tn', 'HCO3tn', 'CO2tv', 'CO2tp'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(co2)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(co2(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('SLATEGRAY')};
            end
        end
    end
    aa = {'MMETt2', 'GLUt2r', 'GLNt2r', 'GLNt7', 'GLUt7', 'LCYSTintv', 'NH4t', 'TCHOLAabcv'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(aa)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(aa(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('ORANGE')};
            end
        end
    end
    cenMet = {'COAtim', 'NAt3_1', 'MALS'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(cenMet)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(cenMet(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('BLUE')};
            end
        end
    end
    lip = {'G3PIt'};
    for i = 1:length(modifiedMap.rxnName)
        a = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        for j = 1:length(lip)
            if ~isempty(a) && string(modifiedMap.rxnName{i}) == string(lip(j))...
                    & string(model2.subSystems(a)) ~= "NonIncluded"
                modifiedMap.rxnColor(i, 1) = {colors('LIME')};
            end
        end
    end
end

if rxnsFluxTask == 5 % change color to dark-color in negative values
    for i = 1:length(modifiedMap.rxnName)
        valueskeys = values(colors)';
        valueskeys (:,2) = keys(colors)';
        a = find(ismember(valueskeys, modifiedMap.rxnColor{i, 1}));
        b = find(ismember(model2.rxns, modifiedMap.rxnName{i}));
        if ~isempty(b) && v1(b) < 0 && string(model2.subSystems(b)) ~= "NonIncluded"
            if string(valueskeys (a,2)) == 'LIME'
                modifiedMap.rxnColor(i, 1) = {colors('GREEN')};
            else
                rxnColor =strcat('DARK', char(valueskeys (a(1,1),2)));
                if string(rxnColor) == 'DARKAQUA'
                    rxnColor = 'DARKCYAN';
                elseif string(rxnColor) == 'DARKFUCHSIA'
                    rxnColor = 'DARKMAGENTA';
                end
                modifiedMap.rxnColor(i, 1) = {colors(rxnColor)};
            end
        end
    end
end

end
