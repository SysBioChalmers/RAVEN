function exportModel(model,fileName,exportToYAML,exportGeneComplexes,supressWarnings)
% exportModel
%   Exports a constraint-based model to an SBML file (L3V1 FBCv2). Optionally
%   writes a human-readable YAML file that can be used for model versioning and
%   comparison only.
%
%   model               a model structure
%   fileName            filename to export the model to (without file extension)
%   exportToYAML        true for additional export of YAML file, false for
%                       only export to SBML file (opt, default false)
%   exportGeneComplexes	true if gene complexes (all gene sets linked with
%                       AND relationship) should be recognised and exported
%                       (opt, default false)
%   supressWarnings     true if warnings should be supressed (opt, default
%                       false)
%
%   The optional YAML file is only meant for comparison between different model
%   versions, e.g. as part of model development on GitHub. SBML is the preferred
%   format for storage and redistribution of models. YAML import is therefore not
%   supported, and the file generated here is not compatible with e.g. COBRApy.
%
%   Usage: exportModel(model,fileName,exportToYAML,exportGeneComplexes,supressWarnings)
%
%   Simonas Marcisauskas, 2018-03-17

if nargin<3
    exportToYAML=false;
end
if nargin<4
    exportGeneComplexes=false;
end
if nargin<5
    supressWarnings=false;
end

% The default SBML format settings, which are used as input for appropriate
% libSBML functions to generate the blank SBML model structure before using
% exporting in with OutputSBML to xml file;
sbmlLevel=3;
sbmlVersion=1;
sbmlPackages = {'fbc','groups'};
sbmlPackageVersions = [2,1];

%Check if the "unconstrained" field is still present. This shows if
%exchange metabolites have been removed
if ~isfield(model,'unconstrained')
    if supressWarnings==false
        EM='There is no unconstrained field in the model structure. This means that no metabolites are considered exchange metabolites';
        dispEM(EM,false);
    end
    model.unconstrained=zeros(numel(model.mets),1);
end

% If model id and name (description) don't exist, make sure that default
% strings are included;
if ~isfield(model,'id')
    fprintf('WARNING: The model is missing the "id" field. Uses "blankID". \n');
    model.id='blankID';
end;
if ~isfield(model,'description')
    fprintf('WARNING: The model is missing the "id" field. Uses "blankName". \n');
    model.description='blankName';
end;

%Check the model structure
if supressWarnings==false
   checkModelStruct(model,false);
end

% Adding several blank fields, if they don't exist already. This is
% to reduce the number of conditions below;
if ~isfield(model,'compMiriams')
    model.compMiriams=cell(numel(model.comps),1);
end;
if ~isfield(model,'inchis')
    model.inchis=cell(numel(model.mets),1);
end;
if ~isfield(model,'metFormulas')
    model.metFormulas=cell(numel(model.mets),1);
end;
if ~isfield(model,'metMiriams')
    model.metMiriams=cell(numel(model.mets),1);
end;

if ~isfield(model,'geneMiriams') && isfield(model,'genes')
    model.geneMiriams=cell(numel(model.genes),1);
end;
if ~isfield(model,'subSystems')
    model.subSystems=cell(numel(model.rxns),1);
end;
if ~isfield(model,'eccodes')
    model.eccodes=cell(numel(model.rxns),1);
end;
if ~isfield(model,'rxnReferences')
    model.rxnReferences=cell(numel(model.rxns),1);
end;
if ~isfield(model,'rxnConfidenceScores')
    model.rxnConfidenceScores=NaN(numel(model.rxns),1);
end;
if ~isfield(model,'rxnNotes')
    model.rxnNotes=cell(numel(model.rxns),1);
end;
if ~isfield(model,'rxnMiriams')
    model.rxnMiriams=cell(numel(model.rxns),1);
end;

if sbmlLevel<3
	%Check if genes have associated compartments
	if ~isfield(model,'geneComps') && isfield(model,'genes')
        if supressWarnings==false
            EM='There are no compartments specified for genes. All genes will be assigned to the first compartment. This is because the SBML structure requires all elements to be assigned to a compartment';
            dispEM(EM,false);
        end
        model.geneComps=ones(numel(model.genes),1);
	end
end

% Convert ids to SBML-convenient format. This is to avoid the data loss
% when unsupported characters are included in ids. Here we are using part
% from convertSBMLID, originating from the COBRA Toolbox;
model.rxns=regexprep(model.rxns,'([^0-9_a-zA-Z])','__${num2str($1+0)}__');
model.mets=regexprep(model.mets,'([^0-9_a-zA-Z])','__${num2str($1+0)}__');
model.comps=regexprep(model.comps,'([^0-9_a-zA-Z])','__${num2str($1+0)}__');
model.genes=regexprep(model.genes,'([^0-9_a-zA-Z])','__${num2str($1+0)}__');

% Generate an empty SBML structure;
modelSBML=getSBMLStructure(sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
modelSBML.metaid=model.id;
modelSBML.id=model.id;
modelSBML.name=model.description;

if isfield(model,'annotation')
    if isfield(model.annotation,'note')
        modelSBML.notes=['<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>',regexprep(model.annotation.note,'<p>|</p>',''),'</p></body></notes>'];        
    end
else
    modelSBML.notes='<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>This file was generated using the exportModel function in RAVEN Toolbox 2.0 and OutputSBML in libSBML </p></body></notes>';
end

modelSBML.annotation=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#meta_' model.id '">'];
if isfield(model,'annotation')
    nameString='';
    if isfield(model.annotation,'familyName')
        if ~isempty(model.annotation.familyName)
            nameString=['<vCard:Family>' model.annotation.familyName '</vCard:Family>'];
        end
    end
    if isfield(model.annotation,'givenName')
        if ~isempty(model.annotation.givenName)
            nameString=[nameString '<vCard:Given>' model.annotation.givenName '</vCard:Given>'];
        end
    end
    email='';
    if isfield(model.annotation,'email')
        if ~isempty(model.annotation.email)
            email=['<vCard:EMAIL>' model.annotation.email '</vCard:EMAIL>'];
        end
    end
    org='';
    if isfield(model.annotation,'organization')
        if ~isempty(model.annotation.organization)
            org=['<vCard:ORG rdf:parseType="Resource"><vCard:Orgname>' model.annotation.organization '</vCard:Orgname></vCard:ORG>'];
        end
    end
    if ~isempty(nameString) || ~isempty(email) || ~isempty(org)
        modelSBML.annotation=[modelSBML.annotation '<dc:creator><rdf:Bag><rdf:li rdf:parseType="Resource">'];
        if ~isempty(nameString)
            modelSBML.annotation=[modelSBML.annotation '<vCard:N rdf:parseType="Resource">' nameString '</vCard:N>'];
        end
        modelSBML.annotation=[modelSBML.annotation email org '</rdf:li></rdf:Bag></dc:creator>'];
    end
end
modelSBML.annotation=[modelSBML.annotation '<dcterms:created rdf:parseType="Resource">'...
'<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'...
'</dcterms:created>'...
'<dcterms:modified rdf:parseType="Resource">'...
'<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'...
'</dcterms:modified>'];

if isfield(model,'annotation')
    if isfield(model.annotation,'taxonomy')
        modelSBML.annotation=[modelSBML.annotation '<bqbiol:is><rdf:Bag><rdf:li rdf:resource="http://identifiers.org/taxonomy/' regexprep(model.annotation.taxonomy,'taxonomy/','') '"/></rdf:Bag></bqbiol:is>'];
    end
end
modelSBML.annotation=[modelSBML.annotation '</rdf:Description></rdf:RDF></annotation>'];

%Prepare compartments
for i=1:numel(model.comps)
    % Adding the default values, as these will be the same in all entries;
    if i==1
        if isfield(modelSBML.compartment, 'sboTerm')
            modelSBML.compartment(i).sboTerm=290;
        end;
        if isfield(modelSBML.compartment, 'spatialDimensions')
            modelSBML.compartment(i).spatialDimensions=3;
        end;
        if isfield(modelSBML.compartment, 'size')
            modelSBML.compartment(i).size=1;
        end;
        if isfield(modelSBML.compartment, 'constant')
            modelSBML.compartment(i).constant=1;
        end;
        if isfield(modelSBML.compartment, 'isSetSize')
            modelSBML.compartment(i).isSetSize=1;
        end;
        if isfield(modelSBML.compartment, 'isSetSpatialDimensions')
            modelSBML.compartment(i).isSetSpatialDimensions=1;
        end;
    end;
    % Copying the default values to the next entry as long as it is not the
    % last one;
    if i<numel(model.comps)
    	modelSBML.compartment(i+1)=modelSBML.compartment(i);
    end;
    
    if isfield(modelSBML.compartment,'metaid')
        modelSBML.compartment(i).metaid=model.comps{i};
    end;
    %Prepare Miriam strings
    if ~isempty(model.compMiriams{i}) && isfield(modelSBML.compartment(i),'annotation')
        modelSBML.compartment(i).annotation=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#meta_' model.comps{i} '">'];
        modelSBML.compartment(i).annotation=[modelSBML.compartment(i).annotation '<bqbiol:is><rdf:Bag>'];
        modelSBML.compartment(i).annotation=[modelSBML.compartment(i).annotation getMiriam(model.compMiriams{i}) '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>']; 
    end;
    if isfield(modelSBML.compartment, 'name')
        modelSBML.compartment(i).name=model.compNames{i};
    end;
    if isfield(modelSBML.compartment, 'id')
        modelSBML.compartment(i).id=model.comps{i};
    end;    

end;

%Begin writing species
for i=1:numel(model.mets)
    % Adding the default values, as these will be the same in all entries;
    if i==1
        if isfield(modelSBML.species, 'sboTerm')
            modelSBML.species(i).sboTerm=247;
        end;
        if isfield(modelSBML.species, 'initialAmount')
            modelSBML.species(i).initialAmount=1;
        end;
        if isfield(modelSBML.species, 'initialConcentration')
            modelSBML.species(i).initialConcentration=0;
        end;
        if isfield(modelSBML.species, 'isSetInitialAmount')
            modelSBML.species(i).isSetInitialAmount=1;
        end;
        if isfield(modelSBML.species, 'isSetInitialConcentration')
            modelSBML.species(i).isSetInitialConcentration=1;
        end;
    end;
    % Copying the default values to the next entry as long as it is not the
    % last one;
    if i<numel(model.mets)
    	modelSBML.species(i+1)=modelSBML.species(i);
    end;
    
    if isfield(modelSBML.species,'metaid')
        modelSBML.species(i).metaid=['M_' model.mets{i}];
    end;
    if isfield(modelSBML.species, 'name')
        modelSBML.species(i).name=model.metNames{i};
    end;
    if isfield(modelSBML.species, 'id')
        modelSBML.species(i).id=['M_' model.mets{i}];
    end;
    if isfield(modelSBML.species, 'compartment')
        modelSBML.species(i).compartment=model.comps{model.metComps(i)};
    end; 
    if isfield(model,'unconstrained')
        if model.unconstrained(i)
            modelSBML.species(i).boundaryCondition=1;
        end;
    end
    if isfield(modelSBML.species, 'fbc_charge') && isfield(model,'metCharges')
        if ~isnan(model.metCharges(i))
            modelSBML.species(i).fbc_charge=model.metCharges(i);
            modelSBML.species(i).isSetfbc_charge=1;
        else
            modelSBML.species(i).isSetfbc_charge=0;
        end;
    end;
    if isfield(modelSBML.species,'annotation')
        if ~isempty(model.metMiriams{i}) || ~isempty(model.metFormulas{i})
            hasInchi=false;
            if ~isempty(model.metFormulas{i})
            %Only export formula if there is no InChI. This is because the
            %metFormulas field is populated by InChIs if available
                if ~isempty(model.inchis{i})
                    hasInchi=true;
                end;
                if hasInchi==false
                    modelSBML.species(i).fbc_chemicalFormula=model.metFormulas{i};
                end
            end
            if ~isempty(model.metMiriams{i}) || hasInchi==true
                modelSBML.species(i).annotation=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#meta_M_' model.mets{i} '">'];
                modelSBML.species(i).annotation=[modelSBML.species(i).annotation '<bqbiol:is><rdf:Bag>'];
                if ~isempty(model.metMiriams{i})
                    modelSBML.species(i).annotation=[modelSBML.species(i).annotation getMiriam(model.metMiriams{i})];
                end;
                if hasInchi==true
                    modelSBML.species(i).annotation=[modelSBML.species(i).annotation '<rdf:li rdf:resource="http://identifiers.org/inchi/InChI=' model.inchis{i} '"/>'];
                end;
                modelSBML.species(i).annotation=[modelSBML.species(i).annotation '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
            end;
        end;
    end;
end; 

if isfield(model,'genes')
    for i=1:numel(model.genes)
        % Adding the default values, as these will be the same in all entries;
        if i==1
            if isfield(modelSBML.fbc_geneProduct, 'sboTerm')
                modelSBML.fbc_geneProduct(i).sboTerm=252;
            end;    
        end;
        % Copying the default values to the next index as long as it is not the
        % last one;
        if i<numel(model.genes)
            modelSBML.fbc_geneProduct(i+1)=modelSBML.fbc_geneProduct(i);
        end;
        
        if isfield(modelSBML.fbc_geneProduct,'metaid')
            modelSBML.fbc_geneProduct(i).metaid=model.genes{i};
        end;
        if ~isempty(model.geneMiriams{i}) && isfield(modelSBML.fbc_geneProduct(i),'annotation')
            modelSBML.fbc_geneProduct(i).annotation=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#meta_' model.genes{i} '">'];
            modelSBML.fbc_geneProduct(i).annotation=[modelSBML.fbc_geneProduct(i).annotation '<bqbiol:is><rdf:Bag>'];
            modelSBML.fbc_geneProduct(i).annotation=[modelSBML.fbc_geneProduct(i).annotation getMiriam(model.geneMiriams{i}) '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>']; 
        end;
        if isfield(modelSBML.fbc_geneProduct, 'fbc_id')
            modelSBML.fbc_geneProduct(i).fbc_id=model.genes{i};
        end;
        if isfield(modelSBML.fbc_geneProduct, 'fbc_label') && isfield(model,'geneShortNames')
            modelSBML.fbc_geneProduct(i).fbc_label=model.geneShortNames{i};
        end;
    end;
    if exportGeneComplexes==true
        %Also add the complexes as genes. This is done by splitting grRules
        %on "or" and adding the ones which contain several genes
        geneComplexes={};
        if isfield(model,'grRules')
            %Only grRules which contain " and " can be complexes
            uniqueRules=unique(model.grRules);
            I=cellfun(@any,strfind(uniqueRules,' and '));
            uniqueRules(~I)=[];
            uniqueRules=strrep(uniqueRules,'(','');
            uniqueRules=strrep(uniqueRules,')','');
            uniqueRules=strrep(uniqueRules,' and ',':');
            for i=1:numel(uniqueRules)
                genes=regexp(uniqueRules(i),' or ','split');
                genes=genes{1}(:);
                %Check which ones are complexes
                I=cellfun(@any,strfind(genes,':'));
                geneComplexes=[geneComplexes;genes(I)];
            end
        end
        geneComplexes=unique(geneComplexes);
        if ~isempty(geneComplexes)
            %Then add them as genes. There is a possiblity that a complex A&B is
            %added as separate from B&A. This isn't really an issue so I don't deal
            %with it
            for i=1:numel(geneComplexes)
                modelSBML.fbc_geneProduct(numel(model.genes)+i)=modelSBML.fbc_geneProduct(1);
                if isfield(modelSBML.fbc_geneProduct,'metaid')
                    modelSBML.fbc_geneProduct(numel(model.genes)+i).metaid=geneComplexes{i};
                end;
                if isfield(modelSBML.fbc_geneProduct,'fbc_id')
                    modelSBML.fbc_geneProduct(numel(model.genes)+i).fbc_id=geneComplexes{i};
                else
                	modelSBML.fbc_geneProduct(i).fbc_label=modelSBML.fbc_geneProduct(i).fbc_id;
                end;           
            end;
        end;
    end;
end;

% Generate a list of unique fbc_bound names
totalValues=[model.lb; model.ub];
totalNames=cell(size(totalValues,1),1);

listUniqueValues=unique(totalValues);

for i=1:length(listUniqueValues)
    listUniqueNames{i,1}=['FB',num2str(i),'N',num2str(abs(round(listUniqueValues(i))))]; % create unique flux bound IDs.
    ind=find(ismember(totalValues,listUniqueValues(i)));
    totalNames(ind)=listUniqueNames(i,1);
end

for i=1:length(listUniqueNames)
    % Adding the default values, as these will be the same in all entries;
	if i==1
        if isfield(modelSBML.parameter, 'constant')
            modelSBML.parameter(i).constant=1;
        end;
        if isfield(modelSBML.parameter, 'isSetValue')
            modelSBML.parameter(i).isSetValue=1;
        end;
	end;
    % Copying the default values to the next index as long as it is not the
    % last one;
    if i<numel(listUniqueNames)
        modelSBML.parameter(i+1)=modelSBML.parameter(i);
    end;
    modelSBML.parameter(i).id=listUniqueNames{i};
    modelSBML.parameter(i).value=listUniqueValues(i);
end

for i=1:numel(model.rxns)
    % Adding the default values, as these will be the same in all entries;
    if i==1
        if isfield(modelSBML.reaction, 'sboTerm')
        	modelSBML.reaction(i).sboTerm=176;
        end;
        if isfield(modelSBML.reaction, 'isSetFast')
        	modelSBML.reaction(i).isSetFast=1;
        end;
    end;
    % Copying the default values to the next index as long as it is not the
    % last one;
    if i<numel(model.rxns)
    	modelSBML.reaction(i+1)=modelSBML.reaction(i);
    end;
    
    if isfield(modelSBML.reaction,'metaid')
    	modelSBML.reaction(i).metaid=['R_' model.rxns{i}];
    end;
    
    % Exporting notes information;
    if (~isnan(model.rxnConfidenceScores(i)) || ~isempty(model.rxnReferences{i}) || ~isempty(model.rxnNotes{i}))
        modelSBML.reaction(i).notes='<notes><body xmlns="http://www.w3.org/1999/xhtml">';
        if ~isnan(model.rxnConfidenceScores(i))
            modelSBML.reaction(i).notes=[modelSBML.reaction(i).notes '<p>Confidence Level: ' num2str(model.rxnConfidenceScores(i)) '</p>'];
        end;
        if ~isempty(model.rxnReferences{i})
            modelSBML.reaction(i).notes=[modelSBML.reaction(i).notes '<p>AUTHORS: ' model.rxnReferences{i} '</p>'];
        end;
        if ~isempty(model.rxnNotes{i})
            modelSBML.reaction(i).notes=[modelSBML.reaction(i).notes '<p>NOTES: ' model.rxnNotes{i} '</p>'];
        end;
        modelSBML.reaction(i).notes=[modelSBML.reaction(i).notes '</body></notes>'];
    end;
    
    % Exporting annotation information from rxnMiriams;
    if (~isempty(model.rxnMiriams{i}) && isfield(modelSBML.reaction(i),'annotation')) || ~isempty(model.eccodes{i})
    	modelSBML.reaction(i).annotation=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#meta_R_' model.rxns{i} '">'];
    	modelSBML.reaction(i).annotation=[modelSBML.reaction(i).annotation '<bqbiol:is><rdf:Bag>'];
        if ~isempty(model.eccodes{i})
        	eccodes=regexp(model.eccodes{i},';','split');
            for j=1:numel(eccodes)
                modelSBML.reaction(i).annotation=[modelSBML.reaction(i).annotation  '<rdf:li rdf:resource="http://identifiers.org/ec-code/' regexprep(eccodes{j},'ec-code/|EC','') '"/>'];
            end;
        end;
    	modelSBML.reaction(i).annotation=[modelSBML.reaction(i).annotation getMiriam(model.rxnMiriams{i}) '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>']; 
    end;
    
    if isfield(modelSBML.reaction, 'name')
        modelSBML.reaction(i).name=model.rxnNames{i};
    end;
    if isfield(modelSBML.reaction, 'id')
        modelSBML.reaction(i).id=['R_' model.rxns{i}];
    end 
    
    % Adding the information about reactants and products;
    involvedMets=addReactantsProducts(model,modelSBML,i);
    for j=1:numel(involvedMets.reactant)
        if j<numel(involvedMets.reactant)
            modelSBML.reaction(i).reactant(j+1)=modelSBML.reaction(i).reactant(j);
        end;
        modelSBML.reaction(i).reactant(j).species=involvedMets.reactant(j).species;
        modelSBML.reaction(i).reactant(j).stoichiometry=involvedMets.reactant(j).stoichiometry;
        modelSBML.reaction(i).reactant(j).isSetStoichiometry=involvedMets.reactant(j).isSetStoichiometry;
        modelSBML.reaction(i).reactant(j).constant=involvedMets.reactant(j).constant;
    end;
    if numel(involvedMets.reactant)==0
        modelSBML.reaction(i).reactant='';
    end;
    for j=1:numel(involvedMets.product)
        if j<numel(involvedMets.product)
            modelSBML.reaction(i).product(j+1)=modelSBML.reaction(i).product(j);
        end;
        modelSBML.reaction(i).product(j).species=involvedMets.product(j).species;
        modelSBML.reaction(i).product(j).stoichiometry=involvedMets.product(j).stoichiometry;
        modelSBML.reaction(i).product(j).isSetStoichiometry=involvedMets.product(j).isSetStoichiometry;
        modelSBML.reaction(i).product(j).constant=involvedMets.product(j).constant;
    end;
    if numel(involvedMets.product)==0
        modelSBML.reaction(i).product='';
    end;
    %Export reversibility information. Reactions are irreversible by
    %default;
    if model.rev(i)==1
        modelSBML.reaction(i).reversible=1;
    end;
    if isfield(model, 'rxnComps')
        modelSBML.reaction(i).compartment=model.comps{model.rxnComps(i)};
    end;
    if isfield(model, 'grRules')
        modelSBML.reaction(i).fbc_geneProductAssociation.fbc_association.fbc_association=model.grRules{i};
    end;
    modelSBML.reaction(i).fbc_lowerFluxBound=totalNames{i};
    modelSBML.reaction(i).fbc_upperFluxBound=totalNames{length(model.lb)+i};
end;

% Prepare subSystems
% Code taken from COBRA functions getModelSubSystems, writeSBML, 
% findRxnsFromSubSystem under
% GNU General Public License v3.0, license file in readme/GPL.MD. Code
% modified for RAVEN.
modelSBML.groups_group.groups_kind = 'partonomy';
modelSBML.groups_group.sboTerm = 633;
tmpStruct=modelSBML.groups_group;
rxns=strcat('R_',model.rxns);
if isfield(model, 'subSystems')
    if ~any(cellfun(@iscell,model.subSystems))
        subSystems = setdiff(model.subSystems,'');
    else
        orderedSubs = cellfun(@(x) columnVector(x),model.subSystems,'UniformOUtput',false);
        subSystems = setdiff(vertcat(orderedSubs{:}),'');
    end
    if isempty(subSystems)
        subSystems = {};
    end
    if ~isempty(subSystems)
    %Build the groups for the group package.
    groupIDs = strcat('group',cellfun(@num2str, num2cell(1:length(subSystems))','UniformOutput',false));    
    for i = 1:length(subSystems)
        cgroup = tmpStruct;
        if ~any(cellfun(@iscell,model.subSystems))
            present = ismember(model.subSystems,subSystems{i});
        else
            present = cellfun(@(x) any(ismember(x,subSystems{i})),model.subSystems);
        end
        groupMembers = rxns(present);
        for j = 1:numel(groupMembers)            
            cMember = tmpStruct.groups_member;
            cMember.groups_idRef = groupMembers{j};
            if j == 1
                cgroup.groups_member = cMember;
            else
                cgroup.groups_member(j) = cMember;
            end
        end
        cgroup.groups_id = groupIDs{i};
        cgroup.groups_name = subSystems{i};
        if i == 1
            modelSBML.groups_group = cgroup;
        else
            modelSBML.groups_group(i) = cgroup;
        end
    end
    end
end
 
% Preparing fbc_objective subfield;

modelSBML.fbc_objective.fbc_type='maximize';
modelSBML.fbc_objective.fbc_id='obj';

ind=find(model.c);

if isempty(ind)
    modelSBML.fbc_objective.fbc_fluxObjective.fbc_coefficient=0;
else
    for i=1:length(ind)
        % Copying the default values to the next index as long as it is not the
        % last one;
        if i<numel(ind)
            modelSBML.reaction(i+1)=modelSBML.reaction(i);
        end;
        values=model.c(model.c~=0);
        modelSBML.fbc_objective(i).fbc_fluxObjective.fbc_reaction=modelSBML.reaction(ind(i)).id;
        modelSBML.fbc_objective(i).fbc_fluxObjective.fbc_coefficient=values(i);
        modelSBML.fbc_objective(i).fbc_fluxObjective.isSetfbc_coefficient=1;
    end;
end;

modelSBML.fbc_activeObjective=modelSBML.fbc_objective.fbc_id;

fbcStr=['http://www.sbml.org/sbml/level', num2str(sbmlLevel), '/version', num2str(sbmlVersion), '/fbc/version',num2str(sbmlPackageVersions(1))];
groupStr=['http://www.sbml.org/sbml/level', num2str(sbmlLevel), '/version', num2str(sbmlVersion), '/groups/version',num2str(sbmlPackageVersions(2))];
modelSBML.namespaces=struct('prefix',{'','fbc','groups'},...
    'uri',{['http://www.sbml.org/sbml/level', num2str(sbmlLevel), '/version', num2str(sbmlVersion), '/core'],...
    fbcStr,groupStr});

if sbmlPackageVersions(1) == 2
    modelSBML.fbc_strict=1;
end;

OutputSBML(modelSBML,strcat(fileName,'.xml'),1,0,[1,0]);

if exportToYAML==true
    % Simplify structure by removing empty fields to reduce YAML size.
    fnm=fieldnames(modelSBML);
    idx=structfun(@isempty,modelSBML);
    modelSBML=rmfield(modelSBML,fnm(idx));
    fnm=fieldnames(modelSBML);
    ws=modelSBML; % Keep working structure, as removing fields below will interfer with indexing
    for i=1:length(fnm)
        if isstruct(ws.(fnm{i}))
            fnm2=fieldnames(ws.(fnm{i}));
            for j=1:length(fnm2)
                if isempty([ws.(fnm{i})(:).(fnm2{j})])
                    modelSBML.(fnm{i})=rmfield(modelSBML.(fnm{i}),fnm2{j});
                else
                    for l=1:length({ws.(fnm{i}).(fnm2{j})})
                        if isstruct(ws.(fnm{i})(l).(fnm2{j}))
                            fnm3=fieldnames(ws.(fnm{i})(l).(fnm2{j}));
                            for k=1:length(fnm3)
                                if isempty([ws.(fnm{i})(l).(fnm2{j})(:).(fnm3{k})])
                                    modelSBML.(fnm{i})(l).(fnm2{j})=rmfield(modelSBML.(fnm{i})(l).(fnm2{j}),fnm3{k});
%                                 else
%                                     isstruct(isstruct({ws.(fnm{i})(l).(fnm2{j}).(fnm3{k})}));
%                                     fnm4=fieldnames(ws.(fnm{i})(l).(fnm2{j}).(fnm3{k}));
%                                     for m=1:length(fnm4)
%                                         if isempty([ws.(fnm{i})(l).(fnm2{j}).(fnm3{k}).(fnm4{m})])
%                                             modelSBML.(fnm{i})(l).(fnm2{j}).(fnm3{k})=rmfield(modelSBML.(fnm{i})(l).(fnm2{j}).(fnm3{k}),fnm4{m});
%                                         end
%                                     end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    YAML.write(strcat(fileName,'.yml'),modelSBML)
end
end


function modelSBML=getSBMLStructure(sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions)
%Returns the blank SBML model structure by using appropriate libSBML
%functions. This creates structure by considering three levels

sbmlFieldNames=getStructureFieldnames('model',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
sbmlDefaultValues=getDefaultValues('model',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);

for i=1:numel(sbmlFieldNames)
	modelSBML.(sbmlFieldNames{1,i})=sbmlDefaultValues{1,i};
    sbmlSubfieldNames=getStructureFieldnames(sbmlFieldNames{1,i},sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
    sbmlSubfieldValues=getDefaultValues(sbmlFieldNames{1,i},sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
    if ~strcmp(sbmlFieldNames{1,i},'event') && ~strcmp(sbmlFieldNames{1,i},'functionDefinition') && ~strcmp(sbmlFieldNames{1,i},'initialAssignment')
        for j=1:numel(sbmlSubfieldNames)
            modelSBML.(sbmlFieldNames{1,i}).(sbmlSubfieldNames{1,j})=sbmlSubfieldValues{1,j};
            sbmlSubsubfieldNames=getStructureFieldnames(sbmlSubfieldNames{1,j},sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
            sbmlSubsubfieldValues=getDefaultValues(sbmlSubfieldNames{1,j},sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
            if ~strcmp(sbmlSubfieldNames{1,j},'modifier') && ~strcmp(sbmlSubfieldNames{1,j},'kineticLaw')
                for k=1:numel(sbmlSubsubfieldNames)
                    % 'compartment' and 'species' field are not supposed to have
                    % their standalone structures if they are subfields or
                    % subsubfields;
                    if ~strcmp(sbmlSubfieldNames{1,j},'compartment') && ~strcmp(sbmlSubfieldNames{1,j},'species')
                        modelSBML.(sbmlFieldNames{1,i}).(sbmlSubfieldNames{1,j}).(sbmlSubsubfieldNames{1,k})=sbmlSubsubfieldValues{1,k};
                    end;   
                    % If it is fbc_association in the third level, we need to
                    % establish the fourth level, since libSBML requires it;
                    if strcmp(sbmlSubsubfieldNames{1,k},'fbc_association')
                        fbc_associationFieldNames=getStructureFieldnames('fbc_association',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
                        fbc_associationFieldValues=getDefaultValues('fbc_association',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
                        for l=1:numel(fbc_associationFieldNames)
                            modelSBML.(sbmlFieldNames{1,i}).(sbmlSubfieldNames{1,j}).(sbmlSubsubfieldNames{1,k}).(fbc_associationFieldNames{1,l})=fbc_associationFieldValues{1,l};
                        end;
                    end;
                end;
            end;
        end;
    end;
    if ~isstruct(modelSBML.(sbmlFieldNames{1,i}))
        modelSBML.(sbmlFieldNames{1,i})=sbmlDefaultValues{1,i};
    end;
end;

modelSBML.unitDefinition.id='mmol_per_gDW_per_hr';

unitFieldNames=getStructureFieldnames('unit',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);
unitDefaultValues=getDefaultValues('unit',sbmlLevel,sbmlVersion,sbmlPackages,sbmlPackageVersions);

kinds={'mole','gram','second'};
exponents=[1 -1 -1];
scales=[-3 0 0];
multipliers=[1 1 1*60*60];

for i=1:numel(unitFieldNames)
    modelSBML.unitDefinition.unit(1).(unitFieldNames{1,i})=unitDefaultValues{1,i};
    for j=1:3
        modelSBML.unitDefinition.unit(j).(unitFieldNames{1,i})=unitDefaultValues{1,i};
        if strcmp(unitFieldNames{1,i},'kind')
            modelSBML.unitDefinition.unit(j).(unitFieldNames{1,i})=kinds{j};
        elseif strcmp(unitFieldNames{1,i},'exponent')
            modelSBML.unitDefinition.unit(j).(unitFieldNames{1,i})=exponents(j);
        elseif strcmp(unitFieldNames{1,i},'scale')
            modelSBML.unitDefinition.unit(j).(unitFieldNames{1,i})=scales(j);
        elseif strcmp(unitFieldNames{1,i},'multiplier')
            modelSBML.unitDefinition.unit(j).(unitFieldNames{1,i})=multipliers(j);
        end;            
    end;
end;
end

function miriamString=getMiriam(miriamStruct)
%Returns a string with list elements for a miriam structure ('<rdf:li
%rdf:resource="http://identifiers.org/go/GO:0005739"/>' for example). This is just
%to speed ut things since this is done many times during the exporting

miriamString='';
if isfield(miriamStruct,'name')
	for i=1:numel(miriamStruct.name)
        miriamString=[miriamString '<rdf:li rdf:resource="http://identifiers.org/' miriamStruct.name{i} '/' miriamStruct.value{i} '"/>'];
	end
end
end

function [tmp_Rxn]=addReactantsProducts(model,sbmlModel,i)
% This function provides reactants and products for particular reaction.
% The function was 'borrowed' from writeSBML in COBRA toolbox, lines
% 663-679;

met_idx = find(model.S(:, i));
tmp_Rxn.product=[];
tmp_Rxn.reactant=[];
for j_met=1:size(met_idx,1)
    tmp_idx = met_idx(j_met,1);
    sbml_tmp_species_ref.species = sbmlModel.species(tmp_idx).id;
    met_stoich = model.S(tmp_idx, i);
    sbml_tmp_species_ref.stoichiometry = abs(met_stoich);
    sbml_tmp_species_ref.isSetStoichiometry=1;
    sbml_tmp_species_ref.constant=1;
    if (met_stoich > 0)
        tmp_Rxn.product = [ tmp_Rxn.product, sbml_tmp_species_ref ];
    else
        tmp_Rxn.reactant = [ tmp_Rxn.reactant, sbml_tmp_species_ref];
    end
end
end

% Code below taken from COBRA Toolbox under GNU General Public License v3.0
% license file in readme/GPL.MD.
function vecT = columnVector(vec)
 
% Converts a vector to a column vector
%
% USAGE:
%
%   vecT = columnVector(vec)
%
% INPUT:
%   vec:     a vector
%
% OUTPUT:
%   vecT:    a column vector
%
% .. Authors:
%     - Original file: Markus Herrgard
%     - Minor changes: Laurent Heirendt January 2017
 
[n, m] = size(vec);

if n < m
    vecT = vec';
else
    vecT = vec;
end
end
 