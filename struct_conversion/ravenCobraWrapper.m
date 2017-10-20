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
%   fields are lost: id, description, annotation, compOutside, compMiriams,
%   rxnComps, geneComps, unconstrained. Boundary metabolites are lost,
%   because COBRA structure does not involve boundary metabolites, so they
%   are removed using simplifyModel before RAVEN -> COBRA conversion. The
%   field 'rev' is also partially lost, but during COBRA -> RAVEN
%   conversion it's reconstructed based on lower bound reaction values
%
%   NOTE: During COBRA -> RAVEN -> COBRA conversion cycle the following
%   fields are lost: b, csense, osense, description, geneEntrezID,
%   metNotes, metSmiles, modelVersion, proteinNames, proteins
%
%   NOTE: The information about mandatory RAVEN fields was taken from
%   checkModelStruct function, whereas the corresponding information about
%   COBRA fields was fetched from verifyModel function
%
%   Usage: newModel=ravenCobraWrapper(model)
%
%   Simonas Marcisauskas, 2017-10-20
%

if isfield(model,'rules')
    isRaven=false;
else
    isRaven=true;
end;

if isRaven
    fprintf('Converting RAVEN structure to COBRA..\n');
    % Converting from RAVEN to COBRA structure;
    
    % Firstly removing boundary metabolites;
    model=simplifyModel(model);

    % Mandatory COBRA fields;
    newModel.rxns=model.rxns;
    newModel.mets=strcat(model.mets,'[',model.comps(model.metComps),']');
    newModel.S=model.S;
    newModel.lb=model.lb;
    newModel.ub=model.ub;
    newModel.c=model.c;
    % b, csense, osense, genes, rules are also mandatory, but defined later
    % to match the order of fields;
    
    % Optional COBRA fields;
    if isfield(model,'rxnNames')
        newModel.rxnNames=model.rxnNames;
    end;
    if isfield(model,'subSystems')
        newModel.subSystems=model.subSystems;
    end;
    if isfield(model,'eccodes')
        newModel.rxnECNumbers=model.eccodes;
    end;
    if isfield(model,'rxnMiriams')
        tmp_rxnkeggid=strrep(extractMiriam(model.rxnMiriams,'kegg.reaction'),'kegg.reaction/','');
        if ~all(cellfun(@isempty,tmp_rxnkeggid))
            newModel.rxnKEGGID=tmp_rxnkeggid;
        end;
    end;
    if isfield(model,'rxnReferences')
        newModel.rxnReferences=model.rxnReferences;
    end;
    if isfield(model,'rxnNotes')    
        newModel.rxnNotes=model.rxnNotes;
    end;
    if isfield(model,'metNames')
        newModel.metNames=model.metNames;
    end;
    if isfield(model,'metFormulas')
        newModel.metFormulas=model.metFormulas;
    end;
    if isfield(model,'metMiriams')
        tmp_kegg_1=strrep(extractMiriam(model.metMiriams,'kegg.compound'),'kegg.compound/','');
        tmp_kegg_2=strrep(extractMiriam(model.metMiriams,'kegg.glycan'),'kegg.glycan/','');
        if ~all(cellfun(@isempty,tmp_kegg_1)) || ~all(cellfun(@isempty,tmp_kegg_2))
            newModel.metKEGGID=regexprep(strcat(tmp_kegg_1, ';',tmp_kegg_2),'^;|;$','');
        end;
        tmp_chebi=strrep(extractMiriam(model.metMiriams,'chebi'),'chebi/','');
        if ~all(cellfun(@isempty,tmp_chebi))
            newModel.metChEBIID=tmp_chebi;
        end;
        tmp_pubchem_1=strrep(extractMiriam(model.metMiriams,'pubchem.compound'),'pubchem.compound/','');
        tmp_pubchem_2=strrep(extractMiriam(model.metMiriams,'pubchem.substance'),'pubchem.substance/','');
        if ~all(cellfun(@isempty,tmp_pubchem_1)) || ~all(cellfun(@isempty,tmp_pubchem_2))
            newModel.metPubChemID=regexprep(strcat(tmp_pubchem_1, ';',tmp_pubchem_2),'^;|;$','');
        end;
    end;
    if isfield(model,'inchis')
    	newModel.metInChIString=regexprep(strcat('InChI=', model.inchis),'^InChI=$','');
    end;
    if isfield(model,'genes')
        newModel.genes=model.genes;
        newModel.rules=grrulesToRules(model);
    else
        fprintf('WARNING: no genes detected. The model therefore may not be exportable to SBML file with writeCbModel\n');
    end;
    if isfield(model,'comps')
        newModel.comps=model.comps;
    end;
    if isfield(model,'compNames')
        newModel.compNames=model.compNames;
    end;
    if isfield(model,'metCharges')
        newModel.metCharges=model.metCharges;
    end;
    if isfield(model,'metMiriams')
        tmp_hmdbid=strrep(extractMiriam(model.metMiriams,'hmdb'),'hmdb/','');
        if ~all(cellfun(@isempty,tmp_hmdbid))
        	newModel.metHMDBID=tmp_hmdbid;
        end;
        tmp_metanetx=strrep(extractMiriam(model.metMiriams,'metanetx.chemical'),'metanetx.chemical/','');
        if ~all(cellfun(@isempty,tmp_metanetx))
            newModel.metMetaNetXID=tmp_metanetx;
        end;
    end;
    newModel.b=zeros(numel(model.mets),1);
    newModel.csense=repmat('E',size(model.mets));
    if isfield(model,'geneMiriams')
        tmp_kegggeneid=strrep(extractMiriam(model.geneMiriams,'kegg.genes'),'kegg.genes/','');
        if ~all(cellfun(@isempty,tmp_kegggeneid))
        	newModel.geneiskegg__46__genesID=tmp_kegggeneid;
        end;
        tmp_genesgdid=strrep(extractMiriam(model.geneMiriams,'sgd'),'sgd/','');
        if ~all(cellfun(@isempty,tmp_genesgdid))
            newModel.geneissgdID=tmp_genesgdid;
        end;
        tmp_proteinuniprotid=strrep(extractMiriam(model.geneMiriams,'uniprot'),'uniprot/','');
        if ~all(cellfun(@isempty,tmp_proteinuniprotid))
            newModel.proteinisuniprotID=tmp_proteinuniprotid;
        end;
    end;
    if isfield(model,'geneShortNames')
        newModel.geneNames=model.geneShortNames;
    end;
    if isfield(model,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=cell2mat(model.rxnConfidenceScores);
    end;
    if isfield(model,'genes')
        newModel.rules=grrulesToRules(model);
    end;
    newModel.osense=-1;
    
    % It seems that grRules, rxnGeneMat and rev are disposable fields in COBRA
    % version, but we export them to make things faster, when converting
    % COBRA structure back to RAVEN;
    if isfield(model,'grRules')
        newModel.grRules=model.grRules;
    end;
    if isfield(model,'rxnGeneMat')
        newModel.rxnGeneMat=model.rxnGeneMat;
    end;
    newModel.rev=model.rev;
 else
    fprintf('Converting COBRA structure to RAVEN..\n');
    % Converting from COBRA to RAVEN structure;
    
    % Mandatory RAVEN fields;
    newModel.rxns=model.rxns;
    newModel.mets=model.mets;
    for i=1:numel(model.comps)
        newModel.mets=regexprep(newModel.mets,['\[', model.comps{i}, '\]$'],'');
        newModel.mets=regexprep(newModel.mets,['\[', model.compNames{i}, '\]$'],'');
    end;
    
    % It some rare cases, there may be overlapping mets due to removal e.g.
    % [c]. To avoid this, we change e.g. [c] into _c;
    if numel(unique(newModel.mets))~=numel(model.mets)
    	newModel.mets=model.mets;
        for i=1:numel(model.comps)
            newModel.mets=regexprep(newModel.mets,'\]$','');
            newModel.mets=regexprep(newModel.mets,['\[', model.comps{i}, '$'],['_', model.comps{i}]);
        end;
    end;
        
    newModel.S=model.S;
    newModel.lb=model.lb;
    newModel.ub=model.ub;
    % Since COBRA no longer contains rev field it is assumed that rxn is
    % reversible if its lower bound is set to zero;
    if ~isfield(model,'rev')
        for i=1:numel(model.rxns)
            if model.lb(i)<0
                newModel.rev(i,1)=1;
            else
                newModel.rev(i,1)=0;
            end;
        end;
    else
    	newModel.rev=model.rev;
    end;
    newModel.c=model.c;
    newModel.b=zeros(numel(model.mets),1);
    if isfield(model,'comps')
        newModel.comps=model.comps;
    else
        % Since 'comps' field is not mandatory in COBRA, it may be required
        % to obtain the non-redundant list of comps from metabolite ids, if
        % 'comps' field is not available;
        newModel.comps=regexprep(model.mets,'^.+\[','');
        newModel.comps=regexprep(newModel.comps,'\]$','');
        newModel.comps=unique(newModel.comps);
    end;
    
    % metComps is also mandatory, but defined later to match the order of fields;
    
    % Fields 'description' and 'id' are also considered as mandatory, but
    % these are added to the model during exportModel/exportToExcelFormat
    % anyway, so there is no point to add this information here;
    
    % Optional RAVEN fields;
    if isfield(model,'compNames')
        newModel.compNames=model.compNames;
    end;
    if isfield(model,'rxnNames')
        newModel.rxnNames=model.rxnNames;
    end;
    if isfield(model,'grRules')
        newModel.grRules=model.grRules;
    else
        model.grRules=rulesTogrrules(model);
        newModel.grRules=model.grRules;
    end;
    if isfield(model,'rxnGeneMat')
        newModel.rxnGeneMat=model.rxnGeneMat;
    elseif isfield(model,'grRules')
        newModel.rxnGeneMat=getRxnGeneMat(model);
    end;
    if isfield(model,'subSystems')
        newModel.subSystems=model.subSystems;
    end;
    if isfield(model,'rxnECNumbers')
        newModel.eccodes=regexprep(model.rxnECNumbers,'EC|EC:','');
    end;
    if isfield(model,'rxnKEGGID')
        for i=1:numel(model.rxns)
            counter=1;
            newModel.rxnMiriams{i,1}=[];
            if ~isempty(model.rxnKEGGID{i})
                newModel.rxnMiriams{i,1}.name{counter,1} = 'kegg.reaction';  
                newModel.rxnMiriams{i,1}.value{counter,1} = model.rxnKEGGID{i};
                counter=counter+1;
            end;
        end;
    end;
    if isfield(model,'rxnNotes')
        newModel.rxnNotes=model.rxnNotes;
    end;
    if isfield(model,'rxnReferences')
        newModel.rxnReferences=model.rxnReferences;
    end;
    if isfield(model,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=num2cell(model.rxnConfidenceScores);
    end;
    if isfield(model,'genes')
        newModel.genes=model.genes;
    end;
    if isfield(model,'geneiskegg__46__genesID') || isfield(model,'geneissgdID') || isfield(model,'metKEGGID')
        for i=1:numel(model.genes)
            counter=1;
            newModel.geneMiriams{i,1}=[];
            if isfield(model,'geneiskegg__46__genesID')
                if ~isempty(model.geneiskegg__46__genesID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'kegg.genes';  
                    newModel.geneMiriams{i,1}.value{counter,1} = model.geneiskegg__46__genesID{i};
                    counter=counter+1;
                end;
            end;
            if isfield(model,'geneissgdID')
                if ~isempty(model.geneissgdID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'sgd';  
                    newModel.geneMiriams{i,1}.value{counter,1} = model.geneissgdID{i};
                    counter=counter+1;
                end;
            end;
            if isfield(model,'proteinisuniprotID')
                if ~isempty(model.proteinisuniprotID{i})
                    newModel.geneMiriams{i,1}.name{counter,1} = 'uniprot';  
                    newModel.geneMiriams{i,1}.value{counter,1} = model.proteinisuniprotID{i};
                    counter=counter+1;
                end;
            end;
        end;
    end;
    if isfield(model,'geneNames')
        newModel.geneShortNames=model.geneNames;
    end;
    newModel.metNames=model.metNames;
    for i=1:numel(model.comps)
        newModel.metNames=regexprep(newModel.metNames,['\[', model.comps{i}, '\]$'],'');
        newModel.metNames=regexprep(newModel.metNames,['\[', model.compNames{i}, '\]$'],'');
    end;
    newModel.metComps=regexprep(model.mets,'^.+\[','');
    newModel.metComps=regexprep(newModel.metComps,'\]$','');
    [~, newModel.metComps]=ismember(newModel.metComps,newModel.comps);
    if isfield(model,'metInChIString')
        newModel.inchis=regexprep(model.metInChIString,'^InChI=','');
    end;
    printWarning=false;
    if isfield(model,'metFormulas')
        newModel.metFormulas=model.metFormulas;
    end;
    if isfield(model,'metChEBIID') || isfield(model,'metHMDBID') || isfield(model,'metKEGGID') || isfield(model,'metPubChemID') || isfield(model,'metMetaNetXID')
        for i=1:numel(model.mets)
            counter=1;
            newModel.metMiriams{i,1}=[];
            if isfield(model,'metChEBIID')
                if ~isempty(model.metChEBIID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'chebi';  
                    newModel.metMiriams{i,1}.value{counter,1} = model.metChEBIID{i};
                    counter=counter+1;
                end;
            end;
            if isfield(model,'metHMDBID')
                if ~isempty(model.metHMDBID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'hmdb';  
                    newModel.metMiriams{i,1}.value{counter,1} = model.metHMDBID{i};
                    counter=counter+1;
                end;
            end;
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
                    end;
                end;
            end;
            if isfield(model,'metPubChemID')
                if ~isempty(model.metPubChemID{i})
                    if strcmp(model.metPubChemID{i}(1:4),'CID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';  
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    elseif strcmp(model.metPubChemID{i}(1:4),'SID:')
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.substance';  
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    else
                        newModel.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';  
                        newModel.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                        printWarning=true;
                    end;
                end;
            end;
            if isfield(model,'metMetaNetXID')
                if ~isempty(model.metMetaNetXID{i})
                    newModel.metMiriams{i,1}.name{counter,1} = 'metanetx.chemical';  
                    newModel.metMiriams{i,1}.value{counter,1} = model.metMetaNetXID{i};
                    counter=counter+1;
                end;
            end;
        end;
    end;
    if printWarning
        fprintf('Could not separate between PubChemID compounds and substances. All annotated PubChemIDs will therefore be assigned as compounds\n');
    end;
    if isfield(model,'metCharges')
        newModel.metCharges=model.metCharges;
    end;
end;

end

function rules=grrulesToRules(model)
    % This function just takes grRules, changes all gene names to
    % 'x(geneNumber)' and also changes 'or' and 'and' relations to
    % corresponding symbols
    replacingGenes=cell([size(model.genes,1) 1]);
    rules=cell([size(model.grRules,1) 1]);
    for i=1:numel(replacingGenes)
        replacingGenes{i}=strcat('x(',num2str(i),')');
    end;
    for i=1:numel(model.grRules)
        rules{i}=regexprep(model.grRules{i},model.genes,replacingGenes);
        rules{i}=regexprep(rules{i},' and ',' & ');
        rules{i}=regexprep(rules{i},' or ',' | ');
    end;
end

function grRules=rulesTogrrules(model)
    % This function just takes rules, changes all gene names to
    % grRules and also changes 'or' and 'and' relations from
    % corresponding symbols
    replacingGenes=cell([size(model.genes,1) 1]);
    grRules=cell([size(model.rules,1) 1]);
    
    for i=1:numel(replacingGenes)
        replacingGenes{i}=strcat('x(',num2str(i),')');
    end;
    for i=1:numel(model.rules)
        grRules{i}=regexprep(model.rules{i},replacingGenes,model.genes);
        grRules{i}=regexprep(grRules{i},' & ',' and ');
        grRules{i}=regexprep(grRules{i},' | ',' or ');
    end;
end

function rxnGeneMat=getRxnGeneMat(model)
    %Check gene association for each reaction and populate rxnGeneMat
    if ~isempty(model.genes)
        rxnGeneMat=zeros(numel(model.rxns),numel(model.genes));
    end
    if ~isempty(model.grRules)
        tempRules=model.grRules;
        for i=1:length(model.rxns)
           %Check that all gene associations have a match in the gene list
           if ~isempty(model.grRules{i})          
               tempRules{i}=regexprep(tempRules{i},' and | or ','>'); %New format: Genes are separated 'and' and 'or' strings with parentheses
               tempRules{i}=regexprep(tempRules{i},'(',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
               tempRules{i}=regexprep(tempRules{i},')',''); %New format: Genes are separated 'and' and 'or' strings with parentheses
               indexesNew=strfind(tempRules{i},'>'); %Old format: Genes are separated by ":" for AND and ";" for OR
               indexes=strfind(tempRules{i},':'); %Old format: Genes are separated by ":" for AND and ";" for OR
               indexes=unique([indexesNew indexes strfind(tempRules{i},';')]);
               if isempty(indexes)
                   %See if you have a match
                   I=find(strcmp(tempRules{i},model.genes));
                   rxnGeneMat(i,I)=1;
               else
                   temp=[0 indexes numel(tempRules{i})+1];
                   for j=1:numel(indexes)+1
                       %The reaction has several associated genes
                       geneName=tempRules{i}(temp(j)+1:temp(j+1)-1);
                       I=find(strcmp(geneName,model.genes));
                       model.rxnGeneMat(i,I)=1;
                   end
               end
           end
        end
    end
    rxnGeneMat=sparse(rxnGeneMat);
end
