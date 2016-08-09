function exportModel(model,fileName,toCOBRA,supressWarnings)
% exportModel
%   Exports a constraint-based model to a SBML file
%
%   model           a model structure
%   fileName        the filename  to export the model to
%   toCOBRA         true if the model should be exported to COBRA format. The
%                   normal format is that which was defined in the Yeast
%                   consensus model (opt, default false)
%   supressWarnings true if warnings should be supressed (opt, default
%                   false)
%
%   Usage: exportModel(model,fileName,toCOBRA,supressWarnings)
%
%   Rasmus Agren, 2013-08-03
%

if nargin<3
    toCOBRA=false;
end
if nargin<4
    supressWarnings=false;
end

%Check if the "unconstrained" field is still present. This shows if
%exchange metabolites have been removed
if ~isfield(model,'unconstrained')
    if supressWarnings==false
        dispEM('There is no unconstrained field in the model structure. This means that no metabolites are considered exchange metabolites',false);
    end
    model.unconstrained=zeros(numel(model.mets),1);
end

%Check the model structure
if supressWarnings==false
   checkModelStruct(model,false); 
end

%For converting illegal characters to their entity reference
model=cleanBadCharsInModel(model);

%Check if genes have associated compartments
if ~isfield(model,'geneComps') && isfield(model,'genes')
    if supressWarnings==false
        dispEM('There are no compartments specified for genes. All genes will be assigned to the first compartment. This is because the SBML structure requires all elements to be assigned to a compartment',false);
    end
    model.geneComps=ones(numel(model.genes),1);
end

%Generate temporary filename
tempF=tempname();

%Open a stream
fid = fopen(tempF,'w');

%Writes the intro
if toCOBRA==false
    intro=['<?xml version="1.0" encoding="UTF-8" ?>'...
    '<sbml xmlns="http://www.sbml.org/sbml/level2/version3" level="2" version="3">'...
    '<model metaid="metaid_' model.id '" id="' model.id '" name="' model.description '">'];
    if isfield(model,'annotation')
        if isfield(model.annotation,'note')
            intro=[intro '<notes><body xmlns="http://www.w3.org/1999/xhtml">' model.annotation.note '</body></notes>\n'];
        end
    else
        intro=[intro '<notes><body xmlns="http://www.w3.org/1999/xhtml">This file was generated using the exportModel function in RAVEN Toolbox</body></notes>\n'];
    end
    intro=[intro '<annotation>'...
    '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
    '<rdf:Description rdf:about="#metaid_' model.id '">'];
    if isfield(model,'annotation')
        nameString='';
        if isfield(model.annotation,'givenName')
            if ~isempty(model.annotation.givenName)
                nameString=['<vCard:Given>' model.annotation.givenName '</vCard:Given>'];
            end
        end
        if isfield(model.annotation,'familyName')
            if ~isempty(model.annotation.familyName)
                nameString=[nameString '<vCard:Family>' model.annotation.familyName '</vCard:Family>'];
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
                org=['<vCard:ORG><vCard:Orgname>' model.annotation.organization '</vCard:Orgname></vCard:ORG>'];
            end
        end
        if ~isempty(nameString) || ~isempty(email) || ~isempty(org)
            intro=[intro '<dc:creator rdf:parseType="Resource"><rdf:Bag><rdf:li rdf:parseType="Resource">'];
            if ~isempty(nameString)
                intro=[intro '<vCard:N rdf:parseType="Resource">' nameString '</vCard:N>'];
            end
            intro=[intro email org '</rdf:li></rdf:Bag></dc:creator>'];
        end
    end
    intro=[intro '<dcterms:created rdf:parseType="Resource">'...
    '<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'... 
    '</dcterms:created>'...
    '<dcterms:modified rdf:parseType="Resource">'...
    '<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'... 
    '</dcterms:modified>'];

    if isfield(model,'annotation')
        if isfield(model.annotation,'taxonomy')
            intro=[intro '<bqbiol:is><rdf:Bag><rdf:li rdf:resource="urn:miriam:' model.annotation.taxonomy '" /></rdf:Bag></bqbiol:is>'];
        end
    end
    intro=[intro '</rdf:Description>'...
    '</rdf:RDF>'...
    '</annotation>'];
else
    intro=['<?xml version="1.0" encoding="UTF-8"?><sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1" xmlns:html="http://www.w3.org/1999/xhtml"><model id="' model.id '" name="' model.description '">'];
    if isfield(model,'annotation')
        if isfield(model.annotation,'note')
            intro=[intro '<notes><body xmlns="http://www.w3.org/1999/xhtml">' model.annotation.note '</body></notes>\n'];
        end
    end
end

intro=[intro '<listOfUnitDefinitions>'...
'<unitDefinition id="mmol_per_gDW_per_hr">'...
'<listOfUnits>'...
'<unit kind="mole" scale="-3"/>'...
'<unit kind="second" multiplier="0.00027778" exponent="-1"/>'...
'</listOfUnits>'...
'</unitDefinition>'...
'</listOfUnitDefinitions>'...
'<listOfCompartments>'];

%Write intro
fprintf(fid,intro);

%Write compartments
for i=1:numel(model.comps)
    %Check if it's outside anything
    if isfield(model, 'compOutside')
        if ~isempty(model.compOutside{i})
            append=[' outside="C_' model.compOutside{i} '" spatialDimensions="3"'];
        else
            append=' spatialDimensions="3"';
        end
    else
        append=' spatialDimensions="3"';
    end
    
    if toCOBRA==false
        fprintf(fid,['<compartment metaid="metaid_C_' model.comps{i} '" id="C_' model.comps{i} '" name="' model.compNames{i} '"' append ' size="1" sboTerm="SBO:0000290">']);

        %Print associated Miriam strings
        if isfield(model,'compMiriams')
            miriamString=getMiriam(model.compMiriams{i});
        else
            miriamString=[];
        end

        if ~isempty(miriamString)
            compinfo=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms='...
                    '"http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
                    'xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#metaid_C_' model.comps{i} '">'...
                    '<bqbiol:is><rdf:Bag>' miriamString '</rdf:Bag></bqbiol:is>'...
                    '</rdf:Description></rdf:RDF></annotation></compartment>'];
        else
            compinfo='</compartment>\n';
        end
        fprintf(fid,compinfo);
    else
        fprintf(fid,['<compartment id="C_' model.comps{i} '" name="' model.compNames{i} '"' append '></compartment>']);
    end
end

intro='</listOfCompartments><listOfSpecies>';
fprintf(fid,intro);

%Begin writing species
for i=1:numel(model.mets)
    if model.unconstrained(i)
        unbounded='true';
    else
        unbounded='false';
    end
    
    if toCOBRA==false
        toprint=['<species metaid="metaid_M_' model.mets{i} '" id="M_' model.mets{i} '" name="' model.metNames{i} '" compartment="C_' model.comps{model.metComps(i)} '" initialAmount="0" boundaryCondition="' unbounded '" sboTerm="SBO:0000299">'];
    else
        %For COBRA format the formula is appended to the metabolite name
        if isfield(model,'metFormulas')
            if ~isempty(model.metFormulas(i))
                append=['_' model.metFormulas{i}];
            else
                append='_';
            end
        else
        	append='_'; 
        end
        toprint=['<species id="M_' model.mets{i} '" name="' model.metNames{i} append '" compartment="C_' model.comps{model.metComps(i)} '" initialAmount="0" boundaryCondition="' unbounded '">'];
    end
    
    %Print some stuff if there is a formula for the compound
    if toCOBRA==false
        if isfield(model,'metFormulas')
            %Only print formula if there is no InChI. This is because the
            %metFormulas field is populated by InChIs if available
            hasInchi=false;
            if isfield(model,'inchis')
                if ~isempty(model.inchis{i})
                    hasInchi=true;
                end
            end

            if ~isempty(model.metFormulas{i}) && hasInchi==false
                toprint=[toprint '<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>FORMULA: '  model.metFormulas{i} '</p></body></notes>'];
            end
        end    
        if isfield(model,'metMiriams')
            miriamString=getMiriam(model.metMiriams{i});
        else
            miriamString='';
        end

        if isfield(model,'inchis')
           if ~isempty(model.inchis{i})
                isInchi=true;
           else
                isInchi=false;
           end
        else
            isInchi=false;
        end

        if any(miriamString) || isInchi==true 
            toprint=[toprint '<annotation>'];

            %Print InChI if available
            if isInchi==true
                toprint=[toprint '<in:inchi xmlns:in="http://biomodels.net/inchi" metaid="metaid_M_' model.mets{i} '_inchi">InChI=' model.inchis{i} '</in:inchi>'];
            end
            %Print some more annotation stuff
            toprint=[toprint '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '...
                'xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
                '<rdf:Description rdf:about="#metaid_M_' model.mets{i} '">'...
                '<bqbiol:is>'...
                '<rdf:Bag>'];
            if isInchi==true
                toprint=[toprint '<rdf:li rdf:resource="#metaid_M_' model.mets{i} '_inchi" />'];
            end

            %Print Miriam and finish up
            toprint=[toprint miriamString '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
        end
    end
    toprint=[toprint '</species>\n'];
    fprintf(fid,toprint);
end

%Write genes and complexes
if toCOBRA==false
    %Write genes
    if isfield(model,'genes')
        for i=1:numel(model.genes)
           toprint=['<species metaid="metaid_E_' int2str(i) '" id="E_' int2str(i) '" name="' model.genes{i} '" compartment="C_' model.comps{model.geneComps(i)} '" initialAmount="0" sboTerm="SBO:0000014">'];

           %Print gene name if present
           if isfield(model,'geneShortNames')
               if ~isempty(model.geneShortNames{i})
                    toprint=[toprint '<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>SHORT NAME: '  model.geneShortNames{i} '</p></body></notes>'];
               end
           end

           if isfield(model,'geneMiriams')
                miriamString=getMiriam(model.geneMiriams{i});
           else
                miriamString=[]; 
           end

           if ~isempty(miriamString)
               toprint=[toprint '<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" '...
                    'xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
                    'xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#metaid_E_' int2str(i) '"><bqbiol:is><rdf:Bag>'...
                    miriamString '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
           end

           toprint=[toprint '</species>\n'];
           fprintf(fid,toprint);
        end

        %Also add the complexes as species. This is done by splitting grRules
        %on "or" and addding the ones which contain several genes
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

        %Then add them as species. There is a possiblity that a complex A&B is
        %added as separate from B&A. This isn't really an issue so I don't deal
        %with it
        for i=1:numel(geneComplexes)
            %The SBO term for the complex is set to the same as for genes. Might not be
            %correct. All complexes are added to the first compartment
            fprintf(fid,['<species metaid="metaid_Cx_' int2str(i) '" id="Cx_' int2str(i) '" name="' geneComplexes{i} '" compartment="C_' model.comps{1} '" initialAmount="0" sboTerm="SBO:0000014"></species>\n']);
        end
    end
end

%Finish metbolites
fprintf(fid,'</listOfSpecies>');

%Add reactions
fprintf(fid,'<listOfReactions>');

if toCOBRA==false
    if isfield(model,'grRules')
        %This is for dealing with complexes
        model.grRules=strrep(model.grRules,'(','');
        model.grRules=strrep(model.grRules,')','');
        model.grRules=strrep(model.grRules,' and ',':');
    end
end

for i=1:length(model.rxns)
    %Get reversibility
    reversible='false';
    if model.rev(i)==1
        reversible='true';
    end
    
    if isfield(model,'subSystems')
        subsystem=model.subSystems{i};
    else
        subsystem=[];
    end
    if isfield(model,'rxnComps')
        rxnComps=model.comps{model.rxnComps(i)};
    else
        rxnComps=[];
    end
    if isfield(model,'eccodes')
        eccode=model.eccodes{i};
    else
        eccode=[];
    end
        
    if toCOBRA==false
        fprintf(fid,['<reaction metaid="metaid_R_' model.rxns{i} '" id="R_' model.rxns{i} '" name="' model.rxnNames{i} '" reversible="' reversible '" sboTerm="SBO:0000176">']);

        if ~isempty(subsystem) || ~isempty(rxnComps)
            toPrint='';
            if ~isempty(subsystem)
               toPrint=[toPrint '<p>SUBSYSTEM: ' subsystem '</p>']; 
            end
            %Compartment isn't an allowed attribute until SBML L3 and I don't
            %feel like changing format just for that
            if ~isempty(rxnComps)
               toPrint=[toPrint '<p>COMPARTMENT: ' rxnComps '</p>']; 
            end
            fprintf(fid,['<notes><body xmlns="http://www.w3.org/1999/xhtml">' toPrint '</body></notes>']);
        end

        %Print annotation
        if isfield(model,'rxnMiriams')
            miriamString=getMiriam(model.rxnMiriams{i});
        else
            miriamString=[];
        end

        %Check to see if there is an associated ec-code
        if any(eccode)
            %Several ec-codes can be delimited by ";"
            eccodes=regexp(eccode,';','split');
            for j=1:numel(eccodes)
                miriamString=[miriamString  '<rdf:li rdf:resource="urn:miriam:ec-code:' eccodes{j} '"/>'];
            end
        end

        if ~isempty(miriamString)
            fprintf(fid,['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '...
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '...
            'xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
            'xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
            '<rdf:Description rdf:about="#metaid_R_' model.rxns{i} '">'...
            '<bqbiol:is>'...
            '<rdf:Bag>'...
            miriamString... 
            '</rdf:Bag>'...
            '</bqbiol:is>'...
            '</rdf:Description>'...
            '</rdf:RDF>'...
            '</annotation>']);
        end
    else
       if isfield(model,'grRules')
            grRules=model.grRules{i};
       else
            grRules=[];
       end
       
       %COBRA format
       fprintf(fid,['<reaction id="R_' model.rxns{i} '" name="' model.rxnNames{i} '" reversible="' reversible '">']);
       
        fprintf(fid,'<notes>');
        if any(grRules)
        	fprintf(fid,['<html:p>GENE_ASSOCIATION: ' grRules '</html:p>']);
        end
        if any(subsystem)
            fprintf(fid,['<html:p>SUBSYSTEM: ' subsystem '</html:p>']);
        end
        if any(eccode)
            fprintf(fid,['<html:p>PROTEIN_CLASS: ' eccode '</html:p>']);
        end
        fprintf(fid,'</notes>');
    end
    
    %The reactants have negative values in the stochiometric matrix
    compounds=model.S(:,i);
    reactants=find(compounds<0);
    products=find(compounds>0);
    
    if any(reactants)
        fprintf(fid,'<listOfReactants>');
        for j=1:length(reactants)
            tempmetname=model.mets{reactants(j)};
            fprintf(fid,['<speciesReference species="M_' tempmetname '" stoichiometry="' num2str(-1*compounds(reactants(j))) '"/>']);
        end  

        fprintf(fid,'</listOfReactants>');
    end
    if any(products)
        fprintf(fid,'<listOfProducts>');
        for j=1:length(products)
           tempmetname=model.mets{products(j)};
           fprintf(fid,['<speciesReference species="M_' tempmetname '" stoichiometry="' num2str(compounds(products(j))) '"></speciesReference>']);
        end  
        fprintf(fid,'</listOfProducts>');
    end
    if toCOBRA==false
        if isfield(model,'genes') && isfield(model,'grRules')
            if ~isempty(model.grRules{i})
                genes=regexp(model.grRules(i),' or ','split');
                genes=genes{1}(:);
                %Check which of them are complexes
                complexes=cellfun(@any,strfind(genes,':'));
                [I normalGenes]=ismember(genes(~complexes),model.genes);
                normalGenes(~I)=[];
                [I complexGenes]=ismember(genes(complexes),geneComplexes);
                complexGenes(~I)=[];
                toPrint='';
                for j=1:numel(normalGenes)
                    toPrint=[toPrint '<modifierSpeciesReference species="E_' num2str(normalGenes(j)) '" />'];
                end
                for j=1:numel(complexGenes)
                    toPrint=[toPrint '<modifierSpeciesReference species="Cx_' num2str(complexGenes(j)) '" />'];
                end
                fprintf(fid,['<listOfModifiers>' toPrint '</listOfModifiers>']);
            end
        end
    end

    %Print constraints, reversibility, and objective. It's assumed that all 
    %information is present in the model structure. This should be ok, since
    %it should have been set to standard values when imported
    
    %Note that the order of parameters is hard-coded. It's the same thing
    %in importModel
    if toCOBRA==false
        fprintf(fid,'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>FLUX_VALUE</ci></math><listOfParameters>');
        fprintf(fid,['<parameter id="LB_R_' model.rxns{i} '" name="LOWER_BOUND" value="' sprintf('%15.8f',model.lb(i)) '" units="mmol_per_gDW_per_hr"/><parameter id="UB_R_' model.rxns{i} '" name="UPPER_BOUND" value="' sprintf('%15.8f',model.ub(i)) '" units="mmol_per_gDW_per_hr"/><parameter id="OBJ_R_' model.rxns{i} '" name="OBJECTIVE_COEFFICIENT" value="' sprintf('%15.8f',model.c(i)) '" units="dimensionless"/><parameter id="FLUX_VALUE" value="0.00000000" units="mmol_per_gDW_per_hr"/>']);
    else
        fprintf(fid,'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><ci> LOWER_BOUND </ci><ci> UPPER_BOUND </ci><ci> OBJECTIVE_COEFFICIENT </ci></apply></math><listOfParameters>');
        fprintf(fid,['<parameter id="LOWER_BOUND" value="' sprintf('%15.8f',model.lb(i)) '"/><parameter id="UPPER_BOUND" value="' sprintf('%15.8f',model.ub(i)) '"/><parameter id="OBJECTIVE_COEFFICIENT" value="' sprintf('%15.8f',model.c(i)) '"/>']);    
    end
    fprintf(fid,'</listOfParameters></kineticLaw>');
    fprintf(fid,'</reaction>\n');
end

%Write outro
outro='</listOfReactions></model></sbml>';
fprintf(fid,outro);

fclose(fid);

%Replace the target file with the temporary file
delete(fileName);
movefile(tempF,fileName);
end

function miriamString=getMiriam(miriamStruct)
%Returns a string with list elements for a miriam structure ('<rdf:li
%rdf:resource="urn:miriam:obo.go:GO:0005739"/>' for example). This is just
%to speed ut things since this is done many times during the exporting

miriamString='';
if isfield(miriamStruct,'name')
	for i=1:numel(miriamStruct.name)
        miriamString=[miriamString '<rdf:li rdf:resource="urn:miriam:' miriamStruct.name{i} ':' miriamStruct.value{i} '"/>'];
	end
end
end
