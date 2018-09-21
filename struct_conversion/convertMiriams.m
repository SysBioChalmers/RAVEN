function model=convertMiriams(model,isRaven);
% convertMiriams
%   Converts MIRIAM annotation between RAVEN and COBRA-style.
%
%   model          RAVEN or COBRA-compatible model structure
%   isRaven        logical whether the model has annotations in RAVEN style
%                  (opt, default checks for presence of met/rxn/geneMiriams
%                  field(s) and uses this as indicator whether model is in
%                  RAVEN style or not)
%
%   model       RAVEN or COBRA-compatible model structure, with
%                  annotation fields in the other style, while
%                  the rest of the model remains in the original style.
%
%   This function bidirectionally converts annotation fields between COBRA
%   and RAVEN style. COBRA has annotations in separate fields (e.g.
%   metKEGGID, metCHEBIID, rxnMetaCycID), while RAVEN has annotations
%   gather in Miriam fields (i.e. rxnMiriams, metMiriams, geneMiriams).
%
%   If isRaven is specified, convertMiriams converts RAVEN -> COBRA if
%   isRaven is true, or COBRA -> RAVEN if isRaven is false. This can be
%   used to force conversion in models whether both types of annotations
%   are present. If no isRaven is specified, convertMiriams checks for
%   presence of met/rxn/geneMiriams field(s) and uses this to set isRaven
%   to true or false.
%
%   Usage: model=convertMiriams(model,isRaven)
%
%   Eduard Kerkhoven, 2018-09-18
%

if nargin<2
    if any(isfield(model,{'rxnMiriams','metMiriams','geneMiriams'}))
        isRaven=true;
    else
        isRaven=false;
    end
end


% Convert towards COBRA style
if isRaven
    if isfield(model,'rxnMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.rxnMiriams);
        miriams=regexprep(miriams,'^[A-Za-z\.]*\/','');
        i=ismember(extractedMiriamNames,'kegg.reaction');
        if any(i)
            model.rxnKEGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'bigg.reaction');
        if any(i)
            model.rxnBIGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'rhea');
        if any(i)
            model.rxnRheaID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metacyc.reaction');
        if any(i)
            model.rxnMetaCycID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'reactome');
        if any(i)
            model.rxnREACTOMEID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sabiork.reaction');
        if any(i)
            model.rxnSABIORKID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'seed.reaction');
        if any(i)
            model.rxnSEEDID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metanetx.reaction');
        if any(i)
            model.rxnMetaNetXID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sbo');
        if any(i)
            model.rxnSBOTerms=miriams(:,i);
        end
        model=rmfield(model,'rxnMiriams');
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
                if i(j) && isfield(model,'metKEGGID')~=1
                    model.metKEGGID=miriams(:,j);
                elseif i(j)
                    model.metKEGGID=strcat(model.metKEGGID,';',miriams(:,j));
                end
            end
            model.metKEGGID=regexprep(model.metKEGGID,'^;|;$','');
        end
        i=ismember(extractedMiriamNames,'chebi');
        if any(i)
            model.metChEBIID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'pubchem');
        if any(i) % Combine Pubchem compounds and substances
            for j=1:length(i)
                if i(j) && isfield(model,'metPubChemID')~=1
                    model.metPubChemID=miriams(:,j);
                elseif i(j)
                    model.metPubChemID=strcat(model.metPubChemID,';',miriams(:,j));
                end
            end
            model.metPubChemID=regexprep(model.metPubChemID,'^;|;$','');
        end
        i=ismember(extractedMiriamNames,'bigg.metabolite');
        if any(i)
            model.metBIGGID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'envipath');
        if any(i)
            model.metEnviPathID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'hmdb');
        if any(i)
            model.metHMDBID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'lipidmaps');
        if any(i)
            model.metLIPIDMAPSID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metacyc.compound');
        if any(i)
            model.metMetaCycID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'reactome.metabolite');
        if any(i)
            model.metREACTOMEID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sabiork.metabolite');
        if any(i)
            model.metSABIORKID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'seed.compound');
        if any(i)
            model.metSEEDID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'slm');
        if any(i)
            model.metSLMID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'metanetx.chemical');
        if any(i)
            model.metMetaNetXID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sbo');
        if any(i)
            model.metSBOTerms=miriams(:,i);
        end
        model=rmfield(model,'metMiriams');
    end
    if isfield(model,'geneMiriams')
        [miriams,extractedMiriamNames]=extractMiriam(model.geneMiriams);
        miriams=regexprep(miriams,'^[A-Za-z\.]*\/','');
        i=ismember(extractedMiriamNames,'kegg.genes');
        if any(i)
            model.geneiskegg__46__genesID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'sgd');
        if any(i)
            model.geneissgdID=miriams(:,i);
        end
        i=ismember(extractedMiriamNames,'uniprot');
        if any(i)
            model.proteinisuniprotID=miriams(:,i);
        end
        model=rmfield(model,'geneMiriams');
    end
else
    % Convert from COBRA style to RAVEN style
    % Define fields
    fields={'rxnBIGGID','rxnKEGGID','rxnMetaCycID','rxnMetaNetXID',...
            'rxnReferences','rxnREACTOMEID','rxnRheaID','rxnSABIORKID',...
            'rxnSBOTerms','rxnSEEDID','geneiskegg__46__genesID',...
            'geneissgdID','proteinisuniprotID','metBIGGID','metChEBIID',...
            'metEnviPathID','metHMDBID',...
            'metKEGGID','metLIPIDMAPSID','metMetaCycID','metMetaNetXID',...
            'metPubChemID','metREACTOMEID','metSABIORKID','metSBOTerms',...
            'metSEEDID','metSLMID'};
    if any(isfield(model,fields(1:10)))
        for i=1:numel(model.rxns)
            counter=1;
            model.rxnMiriams{i,1}=[];
            if isfield(model,'rxnBIGGID')
                if ~isempty(model.rxnBIGGID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'bigg.reaction';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnBIGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnKEGGID')
                if ~isempty(model.rxnKEGGID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'kegg.reaction';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnKEGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnMetaCycID')
                if ~isempty(model.rxnMetaCycID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'metacyc.reaction';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnMetaCycID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnMetaNetXID')
                if ~isempty(model.rxnMetaNetXID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'metanetx.reaction';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnMetaNetXID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnReferences')
                if ~isempty(model.rxnReferences{i})
                    pmids = model.rxnReferences{i};
                    pmids = strsplit(pmids,'; ');
                    for j = 1:length(pmids)
                        model.rxnMiriams{i,1}.name{counter,1} = 'pmid';
                        model.rxnMiriams{i,1}.value{counter,1} = pmids{j};
                        counter=counter+1;
                    end
                end
            end
            if isfield(model,'rxnREACTOMEID')
                if ~isempty(model.rxnREACTOMEID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'reactome';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnREACTOMEID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnRheaID')
                if ~isempty(model.rxnRheaID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'rhea';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnRheaID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSABIORKID')
                if ~isempty(model.rxnSABIORKID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'sabiork.reaction';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnSABIORKID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSBOTerms')
                if ~isempty(model.rxnSBOTerms{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'sbo';
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnSBOTerms{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'rxnSEEDID')
                if ~isempty(model.rxnSEEDID{i})
                    model.rxnMiriams{i,1}.name{counter,1} = 'seed.reaction';
                    %non-official identifiers.org namespace, 'seed'
                    %namespace refers to subsystems
                    model.rxnMiriams{i,1}.value{counter,1} = model.rxnSEEDID{i};
                    counter=counter+1;
                end
            end
        end
    end
    if any(isfield(model,fields(11:13)))
        for i=1:numel(model.genes)
            counter=1;
            model.geneMiriams{i,1}=[];
            if isfield(model,'geneiskegg__46__genesID')
                if ~isempty(model.geneiskegg__46__genesID{i})
                    model.geneMiriams{i,1}.name{counter,1} = 'kegg.genes';
                    model.geneMiriams{i,1}.value{counter,1} = model.geneiskegg__46__genesID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'geneissgdID')
                if ~isempty(model.geneissgdID{i})
                    model.geneMiriams{i,1}.name{counter,1} = 'sgd';
                    model.geneMiriams{i,1}.value{counter,1} = model.geneissgdID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'proteinisuniprotID')
                if ~isempty(model.proteinisuniprotID{i})
                    model.geneMiriams{i,1}.name{counter,1} = 'uniprot';
                    model.geneMiriams{i,1}.value{counter,1} = model.proteinisuniprotID{i};
                    counter=counter+1;
                end
            end
        end
    end
    if any(isfield(model,fields(14:27)))
        for i=1:numel(model.mets)
            counter=1;
            model.metMiriams{i,1}=[];
            if isfield(model,'metBIGGID')
                if ~isempty(model.metBIGGID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'bigg.metabolite';
                    model.metMiriams{i,1}.value{counter,1} = model.metBIGGID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metChEBIID')
                if ~isempty(model.metChEBIID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'chebi';
                    model.metMiriams{i,1}.value{counter,1} = model.metChEBIID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metEnviPathID')
                if ~isempty(model.metChEBIID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'envipath';%not in identifiers.org
                    model.metMiriams{i,1}.value{counter,1} = model.metEnviPathID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metHMDBID')
                if ~isempty(model.metHMDBID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'hmdb';
                    model.metMiriams{i,1}.value{counter,1} = model.metHMDBID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metKEGGID')
                if ~isempty(model.metKEGGID{i})
                    if strcmp(model.metKEGGID{i}(1),'C')
                        model.metMiriams{i,1}.name{counter,1} = 'kegg.compound';
                        model.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    elseif strcmp(model.metKEGGID{i}(1),'G')
                        model.metMiriams{i,1}.name{counter,1} = 'kegg.glycan';
                        model.metMiriams{i,1}.value{counter,1} = model.metKEGGID{i};
                        counter=counter+1;
                    end
                end
            end
            if isfield(model,'metLIPIDMAPSID')
                if ~isempty(model.metLIPIDMAPSID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'lipidmaps';
                    model.metMiriams{i,1}.value{counter,1} = model.metLIPIDMAPSID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metMetaCycID')
                if ~isempty(model.metMetaCycID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'metacyc.compound';
                    model.metMiriams{i,1}.value{counter,1} = model.metMetaCycID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metMetaNetXID')
                if ~isempty(model.metMetaNetXID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'metanetx.chemical';
                    model.metMiriams{i,1}.value{counter,1} = model.metMetaNetXID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metPubChemID')
                if ~isempty(model.metPubChemID{i})
                    if length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'CID:')
                        model.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        model.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    elseif length(model.metPubChemID{i})>3 && strcmp(model.metPubChemID{i}(1:4),'SID:')
                        model.metMiriams{i,1}.name{counter,1} = 'pubchem.substance';
                        model.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                    else
                        model.metMiriams{i,1}.name{counter,1} = 'pubchem.compound';
                        model.metMiriams{i,1}.value{counter,1} = model.metPubChemID{i};
                        counter=counter+1;
                        printWarning=true;
                    end
                end
            end
            if isfield(model,'metREACTOMEID')
                if ~isempty(model.metREACTOMEID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'reactome.metabolite';
                    %non-official identifiers.org, 'reactome' namespace
                    %refers to reactions
                    model.metMiriams{i,1}.value{counter,1} = model.metREACTOMEID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSABIORKID')
                if ~isempty(model.metSEEDID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'sabiork.metabolite';
                    %non-official identifiers.org namespace,
                    %'sabiork.reaction' refers to reactions
                    model.metMiriams{i,1}.value{counter,1} = model.metSABIORKID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSBOTerms')
                if ~isempty(model.metSBOTerms{i})
                    model.metMiriams{i,1}.name{counter,1} = 'sbo';
                    model.metMiriams{i,1}.value{counter,1} = model.metSBOTerms{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSEEDID')
                if ~isempty(model.metSEEDID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'seed.compound';
                    model.metMiriams{i,1}.value{counter,1} = model.metSEEDID{i};
                    counter=counter+1;
                end
            end
            if isfield(model,'metSLMID')
                if ~isempty(model.metSLMID{i})
                    model.metMiriams{i,1}.name{counter,1} = 'swisslipid';
                    model.metMiriams{i,1}.value{counter,1} = model.metSLMID{i};
                    counter=counter+1;
                end
            end
        end
    end
    if exist('printWarning') && printWarning
        fprintf('Could not determine whether PubChemIDs are compounds (CID)\n or substances (SID). All annotated PubChemIDs will therefore \n be assigned as compounds (CID).\n');
    end
    oldFields=isfield(model,fields);
    model=rmfield(model,fields(oldFields));        
end


