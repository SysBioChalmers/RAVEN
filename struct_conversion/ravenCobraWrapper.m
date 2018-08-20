function newModel=ravenCobraWrapper(model)
% ravenCobraWrapper
%   Converts between RAVEN and COBRA structures
%
%   model          a RAVEN/COBRA-compatible model structure
%
%   newModel       a COBRA/RAVEN-compatible model structure
%   
%   This function is a bidirectional tool to convert between RAVEN and COBRA
%   structures. It recognises COBRA structure by checking field 'rules'
%   existense, which is only found in COBRA Toolbox structure.
%
%   NOTE: During RAVEN -> COBRA -> RAVEN conversion cycle the following
%   fields are lost: annotation, compOutside, compMiriams,
%   rxnComps, geneComps, unconstrained. Boundary metabolites are lost,
%   because COBRA structure does not involve boundary metabolites, so they
%   are removed using simplifyModel before RAVEN -> COBRA conversion. The
%   field 'rev' is also partially lost, but during COBRA -> RAVEN
%   conversion it's reconstructed based on lower bound reaction values
%
%   NOTE: During COBRA -> RAVEN -> COBRA conversion cycle the following
%   fields are lost: b, csense, osenseStr, description, geneEntrezID,
%   metNotes, metSmiles, modelVersion, proteinNames, proteins
%
%   NOTE: The information about mandatory RAVEN fields was taken from
%   checkModelStruct function, whereas the corresponding information about
%   COBRA fields was fetched from verifyModel function
%
%   Usage: newModel=ravenCobraWrapper(model)
%
%   Benjamín J. Sánchez, 2018-08-13
%

if isfield(model,'rules')
    isRaven=false;
else
    isRaven=true;
end

if isRaven
    fprintf('Converting RAVEN structure to COBRA..\n');
    %Convert from RAVEN to COBRA structure
    
    %Firstly remove boundary metabolites
    model=simplifyModel(model);
    
    %Mandatory COBRA fields
    newModel.rxns=model.rxns;
    newModel.mets=strcat(model.mets,'[',model.comps(model.metComps),']');
    newModel.S=model.S;
    newModel.lb=model.lb;
    newModel.ub=model.ub;
    newModel.c=model.c;
    %b, csense, osenseStr, genes, rules are also mandatory, but defined
    %later to match the order of fields
    
    %Optional COBRA fields
    if isfield(model,'id')
        newModel.modelID=model.id;
    end
    if isfield(model,'description')
        newModel.modelName=model.description;
    end
    if isfield(model,'rxnNames')
        newModel.rxnNames=model.rxnNames;
    end
    if isfield(model,'subSystems')
        newModel.subSystems=model.subSystems;
    end
    if isfield(model,'eccodes')
        newModel.rxnECNumbers=model.eccodes;
    end
    if isfield(model,'rxnMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.rxnMiriams);
        miriams=regexprep(miriams,'^[A-Za-z\.]*\/','');
        i=ismember(extractedMiriamNames,'kegg.reaction');
        if any(i)
            newModel.rxnKEGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'bigg.reaction');
        if any(i)
            newModel.rxnBIGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'rhea');
        if any(i)
            newModel.rxnRheaID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metacyc.reaction');
        if any(i)
            newModel.rxnMetaCycID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'reactome');
        if any(i)
            newModel.rxnREACTOMEID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sabiork.reaction');
        if any(i)
            newModel.rxnSABIORKID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'seed.reaction');
        if any(i)
            newModel.rxnSEEDID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metanetx.reaction');
        if any(i)
            newModel.rxnMetaNetXID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sbo');
        if any(i)
            newModel.rxnSBOTerms=miriams(:,i);
        end
    end
    if isfield(model,'rxnReferences')
        newModel.rxnReferences=model.rxnReferences;
    end
    if isfield(model,'rxnNotes')
        newModel.rxnNotes=model.rxnNotes;
    end
    if isfield(model,'metNames')
        newModel.metNames=strcat(model.metNames,' [',model.compNames(model.metComps),']');
    end
    if isfield(model,'metFormulas')
        newModel.metFormulas=model.metFormulas;
    end
    if isfield(model,'metMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.metMiriams);
        miriams=regexprep(miriams,'^[A-Za-z\.]*\/','');
        %Shorten miriam names for KEGG and PubChem. These shorter names
        %will be used later to concatenate KEGG COMPOUND/GLYCAN and PubChem
        %Compound/Substance, into corresponding COBRA model fields
        extractedMiriamNames=regexprep(extractedMiriamNames,'^kegg\..+','kegg');
        extractedMiriamNames=regexprep(extractedMiriamNames,'^pubchem\..+','pubchem');
        i=ismember(extractedMiriamNames,'kegg');
        if any(i) % Combine KEGG compounds and glycans
            for j=1:length(i)
                if i(j) && isfield(newModel,'metKEGGID')~=1
                    newModel.metKEGGID=miriams(:,j);
                elseif i(j)
                    newModel.metKEGGID=strcat(newModel.metKEGGID,';',miriams(:,j));
                end
            end
            newModel.metKEGGID=regexprep(newModel.metKEGGID,'^;|;$','');
        end
        i=ismember(extractedMiriamNames,'chebi');
        if any(i)
            newModel.metChEBIID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'pubchem');
        if any(i) % Combine Pubchem compounds and substances
            for j=1:length(i)
                if i(j) && isfield(newModel,'metPubChemID')~=1
                    newModel.metPubChemID=miriams(:,j);
                elseif i(j)
                    newModel.metPubChemID=strcat(newModel.metPubChemID,';',miriams(:,j));
                end
            end
            newModel.metPubChemID=regexprep(newModel.metPubChemID,'^;|;$','');
        end
        i=ismember(extractedMiriamNames,'bigg.metabolite');
        if any(i)
            newModel.metBIGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'envipath');
        if any(i)
            newModel.metEnviPathID=miriams(:,i);
        end        
        i=ismember(extractedMiriamNames,'hmdb');
        if any(i)
            newModel.metHMDBID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'lipidmaps');
        if any(i)
            newModel.metLIPIDMAPSID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metacyc.compound');
        if any(i)
            newModel.metMetaCycID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'reactome.metabolite');
        if any(i)
            newModel.metREACTOMEID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sabiork.metabolite');
        if any(i)
            newModel.metSABIORKID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'seed.compound');
        if any(i)
            newModel.metSEEDID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'slm');
        if any(i)
            newModel.metSLMID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metanetx.chemical');
        if any(i)
            newModel.metMetaNetXID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sbo');
        if any(i)
            newModel.metSBOTerms=miriams(:,i);
        end
    end
    if isfield(model,'inchis')
        newModel.metInChIString=regexprep(strcat('InChI=', model.inchis),'^InChI=$','');
    end
    if isfield(model,'genes')
        newModel.genes=model.genes;
        newModel.rules=grrulesToRules(model);
    else
        fprintf('WARNING: no genes detected. The model therefore may not be exportable to SBML file with writeCbModel\n');
    end
    if isfield(model,'comps')
        newModel.comps=model.comps;
    end
    if isfield(model,'compNames')
        newModel.compNames=model.compNames;
    end
    if isfield(model,'metCharges')
        newModel.metCharges=model.metCharges;
    end
    newModel.b=zeros(numel(model.mets),1);
    newModel.csense=repmat('E',size(model.mets));
    if isfield(model,'geneMiriams')
       [miriams,extractedMiriamNames]=extractMiriam(model.geneMiriams);
        miriams=regexprep(miriams,'^[A-Za-z\.]*\/','');
        i=ismember(extractedMiriamNames,'kegg.genes');
        if any(i)
            newModel.geneiskegg__46__genesID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'kegg.genes');
        if any(i)
            newModel.geneiskegg__46__genesID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sgd');
        if any(i)
            newModel.geneissgdID=miriams(:,i);
        end      
        i=ismember(extractedMiriamNames,'uniprot');
        if any(i)
            newModel.proteinisuniprotID=miriams(:,i);
        end      
    end
    if isfield(model,'geneShortNames')
        newModel.geneNames=model.geneShortNames;
    end
    if isfield(model,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=model.rxnConfidenceScores;
    end
    if isfield(model,'genes')
        newModel.rules=grrulesToRules(model);
    end
    newModel.osenseStr='max';
    
    %It seems that grRules, rxnGeneMat and rev are disposable fields in
    %COBRA version, but we export them to make things faster, when
    %converting COBRA structure back to RAVEN
    if isfield(model,'rxnGeneMat')
        newModel.rxnGeneMat=model.rxnGeneMat;
    end
    if isfield(model,'grRules')
        [grRules, rxnGeneMat] = standardizeGrRules(model,true);
        newModel.grRules      = grRules;
        %Incorporate a rxnGeneMat consistent with standardized grRules
        newModel.rxnGeneMat   = rxnGeneMat;
    end
    newModel.rev=model.rev;
else
    fprintf('Converting COBRA structure to RAVEN..\n');
    %Convert from COBRA to RAVEN structure
    
    %Mandatory RAVEN fields
    newModel.rxns=model.rxns;
    newModel.mets=model.mets;
    for i=1:numel(model.comps)
        newModel.mets=regexprep(newModel.mets,['\[', model.comps{i}, '\]$'],'');
        newModel.mets=regexprep(newModel.mets,['\[', model.compNames{i}, '\]$'],'');
    end
    
    %In some rare cases, there may be overlapping mets due to removal e.g.
    %[c]. To avoid this, we change e.g. [c] into _c
    if numel(unique(newModel.mets))~=numel(model.mets)
        newModel.mets=model.mets;
        for i=1:numel(model.comps)
            newModel.mets=regexprep(newModel.mets,'\]$','');
            newModel.mets=regexprep(newModel.mets,['\[', model.comps{i}, '$'],['_', model.comps{i}]);
        end
    end
    
    newModel.S=model.S;
    newModel.lb=model.lb;
    newModel.ub=model.ub;
    %Since COBRA no longer contains rev field it is assumed that rxn is
    %reversible if its lower bound is set to zero
    if ~isfield(model,'rev')
        for i=1:numel(model.rxns)
            if model.lb(i)<0
                newModel.rev(i,1)=1;
            else
                newModel.rev(i,1)=0;
            end
        end
    else
        newModel.rev=model.rev;
    end
    newModel.c=model.c;
    newModel.b=zeros(numel(model.mets),1);
    if isfield(model,'comps')
        newModel.comps=model.comps;
    else
        %Since 'comps' field is not mandatory in COBRA, it may be required
        %to obtain the non-redundant list of comps from metabolite ids, if
        %'comps' field is not available
        newModel.comps=regexprep(model.mets,'^.+\[','');
        newModel.comps=regexprep(newModel.comps,'\]$','');
        newModel.comps=unique(newModel.comps);
    end
    
    %metComps is also mandatory, but defined later to match the order of
    %fields
    
    %Fields 'description' and 'id' are also considered as mandatory, but
    %these are added to the model during exportModel/exportToExcelFormat
    %anyway, so there is no point to add this information here
    
    %Optional RAVEN fields
    if isfield(model,'modelID')
        newModel.id=model.modelID;
    end
    if isfield(model,'modelName')
        newModel.description=model.modelName;
    end
    if isfield(model,'compNames')
        newModel.compNames=model.compNames;
    end
    if isfield(model,'rxnNames')
        newModel.rxnNames=model.rxnNames;
    end
    if isfield(model,'grRules')
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        newModel.grRules     = grRules;
        newModel.rxnGeneMat  = rxnGeneMat;
    else
        model.grRules        = rulesTogrrules(model);
        [grRules,rxnGeneMat] = standardizeGrRules(model,true);
        newModel.grRules     = grRules;
        newModel.rxnGeneMat  = rxnGeneMat;
    end
    if isfield(model,'subSystems')
        newModel.subSystems=model.subSystems;
    end
    if isfield(model,'rxnECNumbers')
        newModel.eccodes=regexprep(model.rxnECNumbers,'EC|EC:','');
    end
    if any(isfield(model,{'rxnKEGGID','rxnReferences','rxnBIGGID',...
            'rxnMetaCycID','rxnREACTOMEID','rxnRheaID','rxnSABIORKID',...
            'rxnSEEDID','rxnMetaNetXID','rxnSBOTerms'}))
        for i=1:numel(model.rxns)
            counter=1;
            newModel.rxnMiriams{i,1}=[];
            if isfield(model,'rxnKEGGID')
                if ~isempty(model.rxnKEGGID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'kegg.reaction';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnKEGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnBIGGID')
                if ~isempty(model.rxnBIGGID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'bigg.reaction';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnBIGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnMetaCycID')
                if ~isempty(model.rxnMetaCycID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'metacyc.reaction';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnMetaCycID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnREACTOMEID')
                if ~isempty(model.rxnREACTOMEID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'reactome';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnREACTOMEID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnRheaID')
                if ~isempty(model.rxnRheaID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'rhea';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnRheaID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSABIORKID')
                if ~isempty(model.rxnSABIORKID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'sabiork.reaction';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnSABIORKID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSEEDID')
                if ~isempty(model.rxnSEEDID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'seed.reaction';
                    %non-official identifiers.org namespace, 'seed'
                    %namespace refers to subsystems
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnSEEDID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnMetaNetXID')
                if ~isempty(model.rxnMNXID{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'metanetx.reaction';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnMetaNetXID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSBOTerms')
                if ~isempty(model.rxnSBOTerms{i})
                    newModel.rxnMiriams{i,1}.name{counter,1} = 'sbo';
                    newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnSBOTerms{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnReferences')
                if ~isempty(model.rxnReferences{i})
                    pmids = model.rxnReferences{i};
                    pmids = strsplit(pmids,'; ');
                    for j = 1:length(pmids)
                        newModel.rxnMiriams{i,1}.name{counter,1} = 'pmid';
                        newModel.rxnMiriams{i,1}.value{counter,1} = pmids{j};
                        counter=counter+1;
                    end
                end
            end
        end
    end
    if isfield(model,'rxnReferences')
        %if a rxnReferences field is all numeric, it's a pubmedID, and then
        %it's already in rxnMiriams
        %implement regexp if we want to filter those out
        newModel.rxnReferences=model.rxnReferences;
    end
    if isfield(model,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=model.rxnConfidenceScores;
    end
    if isfield(model,'genes')
        newModel.genes=model.genes;
    end
    if isfield(model,'geneiskegg__46__genesID') || isfield(model,'geneissgdID') || isfield(model,'metKEGGID')
        for i=1:numel(model.genes)
            counter=1;
            newModel.geneMiriams{i,1}=[];
            if isfield(model,'geneiskegg__46__genesID')
                if ~isempty(model.geneiskegg__46__genesID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'kegg.genes';
                    newModel.geneMiriams{i,1}.value{counter,1} = model.geneiskegg__46__genesID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'geneissgdID')
                if ~isempty(model.geneissgdID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'sgd';
                    newModel.geneMiriams{i,1}.value{counter,1} = model.geneissgdID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'proteinisuniprotID')
                if ~isempty(model.proteinisuniprotID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'uniprot';
                    newModel.geneMiriams{i,1}.value{counter,1} = model.proteinisuniprotID{i};
                    counter=counter+1;
                end
            end
        end
    end
    if isfield(model,'geneNames')
        newModel.geneShortNames=model.geneNames;
    end
    newModel.metNames=model.metNames;
    for i=1:numel(model.comps)
        newModel.metNames=regexprep(newModel.metNames,['\[', model.comps{i}, '\]$'],'');
        newModel.metNames=regexprep(newModel.metNames,['\[', model.compNames{i}, '\]$'],'');
    end
    newModel.metNames=deblank(newModel.metNames);
    newModel.metComps=regexprep(model.mets,'^.+\[','');
    newModel.metComps=regexprep(newModel.metComps,'\]$','');
    [~, newModel.metComps]=ismember(newModel.metComps,newModel.comps);
    if isfield(model,'metInChIString')
        newModel.inchis=regexprep(model.metInChIString,'^InChI=','');
    end
    printWarning=false;
    if isfield(model,'metFormulas')
        newModel.metFormulas=model.metFormulas;
    end
    if any(isfield(model,{'metChEBIID','metEnviPathID','metHMDBID','metKEGGID',...
            'metPubChemID','metMetaNetXID','metBIGGID','metLIPIDMAPSID',...
            'metMetaCycID','metREACTOMEID','metSABIORKID','metSEEDID',...
            'metSLMID','metSBOTerms'}))
        for i=1:numel(model.mets)
            counter=1;
            newModel.metMiriams{i,1}=[];
            if isfield(model,'metChEBIID')
                if ~isempty(model.metChEBIID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'chebi';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metChEBIID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metEnviPathID')
                if ~isempty(model.metChEBIID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'envipath';%not in identifiers.org
                    newModel.metMiriams{i,1}.value{counter,1} = model.metEnviPathID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metHMDBID')
                if ~isempty(model.metHMDBID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'hmdb';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metHMDBID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metKEGGID')
                if ~isempty(model.metKEGGID{i})
                    if strcmp(model.metKEGGID{i}(1),'C')
                        newModel.metMiriams{i,1}.name{counter,1} = 'kegg.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    elseif strcmp(model.metKEGGID{i}(1),'G')
                        newModel.metMiriams{i,1}.name{counter,1} = 'kegg.glycan';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    end
                end
            end
            if isfield(model,'metPubChemID')
                if ~isempty(model.metPubChemID{i})
                    if length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'CID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    elseif length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'SID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.substance';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    else
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                        printWarning=true;
                    end
                end
            end
            if isfield(model,'metMetaNetXID')
                if ~isempty(model.metMetaNetXID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'metanetx.chemical';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metMetaNetXID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metBiGGID')
                if ~isempty(model.metBiGGID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'bigg.metabolite';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metBiGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metLIPIDMAPSID')
                if ~isempty(model.metLIPIDMAPSID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'lipidmaps';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metLIPIDMAPSID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metMetaCycID')
                if ~isempty(model.metMetaCycID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'metacyc.compound';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metMetaCycID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metREACTOMEID')
                if ~isempty(model.metREACTOMEID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'reactome.metabolite';
                    %non-official identifiers.org, 'reactome' namespace
                    %refers to reactions
                    newModel.metMiriams{i,1}.value{counter,1} = model.metREACTOMEID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSEEDID')
                if ~isempty(model.metSEEDID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'seed.compound';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metSEEDID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSABIORKID')
                if ~isempty(model.metSEEDID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'sabiork.metabolite';
                    %non-official identifiers.org namespace,
                    %'sabiork.reaction' refers to reactions
                    newModel.metMiriams{i,1}.value{counter,1} = model.metSABIORKID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSLMID')
                if ~isempty(model.metSLMID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'swisslipid';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metSLMID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSBOTerms')
                if ~isempty(model.metSBOTerms{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'sbo';
                    newModel.metMiriams{i,1}.value{counter,1} = model.metSBOTerms{i};
                    counter=counter+1;
                end
            end
        end
    end
    if printWarning
        fprintf('Could not determine whether PubChemIDs are compounds (CID)\n or substances (SID). All annotated PubChemIDs will therefore \n be assigned as compounds (CID).\n');
    end
    if isfield(model,'metCharges')
        newModel.metCharges=model.metCharges;
    end
end
end

function rules=grrulesToRules(model)
%This function just takes grRules, changes all gene names to
%'x(geneNumber)' and also changes 'or' and 'and' relations to corresponding
%symbols
replacingGenes=cell([size(model.genes,1) 1]);
for i=1:numel(replacingGenes)
    replacingGenes{i}=strcat('x(',num2str(i),')');
end
rules = strcat({' '},model.grRules,{' '});
for i=1:length(model.genes)
    rules=regexprep(rules,[' ' model.genes{i} ' '],[' ' replacingGenes{i} ' ']);
    rules=regexprep(rules,['(' model.genes{i} ' '],['(' replacingGenes{i} ' ']);
    rules=regexprep(rules,[' ' model.genes{i} ')'],[' ' replacingGenes{i} ')']);
end
rules=regexprep(rules,' and ',' & ');
rules=regexprep(rules,' or ',' | ');
rules=strtrim(rules);
end

function grRules=rulesTogrrules(model)
%This function takes rules, replaces &/| for and/or, replaces the x(i)
%format with the actual gene ID, and takes out extra whitespace and
%redundant parenthesis introduced by COBRA, to create grRules.
grRules = strrep(model.rules,'&','and');
grRules = strrep(grRules,'|','or');
for i = 1:length(model.genes)
    grRules = strrep(grRules,['x(' num2str(i) ')'],model.genes{i});
end
grRules = strrep(grRules,'( ','(');
grRules = strrep(grRules,' )',')');
grRules = regexprep(grRules,'^(',''); %rules that start with a "("
grRules = regexprep(grRules,')$',''); %rules that end with a ")"
end
