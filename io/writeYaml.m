function writeYaml(model,name)
% writeYaml
%   Writes a yaml file matching (roughly) the cobrapy yaml structure
%
%   model       a model structure
%   name        name that the file will have
%
%   Usage: writeYaml(model,name)
%
%   Benjamin Sanchez, 2018-04-12
%

%Check that model is in RAVEN format:
if isfield(model,'rules')
    model = ravenCobraWrapper(model);
end

%Simplify Miriam fields:
if isfield(model,'metMiriams')
    model.metCHEBI = extractMiriam(model.metMiriams,'chebi');
    model.metCHEBI = strrep(model.metCHEBI,'chebi/','');
    model.metKEGG  = extractMiriam(model.metMiriams,'kegg.compound');
    model.metKEGG  = strrep(model.metKEGG,'kegg.compound/','');
end
if isfield(model,'rxnMiriams')
    model.rxnKEGG  = extractMiriam(model.rxnMiriams,'kegg.reaction');
    model.rxnKEGG  = strrep(model.rxnKEGG,'kegg.reaction/','');
end

%Open file:
fid = fopen(name,'wt');
fprintf(fid,'!!omap\n');

%Metabolites:
fprintf(fid,'- metabolites:\n');
[~,pos] = sort(model.mets);
for i = 1:length(model.mets)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'mets',        'txt', pos(i), '- id')
    writeField(model, fid, 'metNames',    'txt', pos(i), '- name')
    writeField(model, fid, 'metComps',    'txt', pos(i), '- compartment')
    writeField(model, fid, 'metFormulas', 'txt', pos(i), '- formula')
    writeField(model, fid, 'metCharges',  'num', pos(i), '- charge')    
    writeField(model, fid, 'metMiriams',  'txt', pos(i), '- annotation')
end

%Reactions:
fprintf(fid,'- reactions:\n');
[~,pos] = sort(model.rxns);
for i = 1:length(model.rxns)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'rxns',                'txt', pos(i), '- id')
    writeField(model, fid, 'rxnNames',            'txt', pos(i), '- name')
    writeField(model, fid, 'S',                   'txt', pos(i), '- metabolites')
    writeField(model, fid, 'lb',                  'num', pos(i), '- lower_bound')
    writeField(model, fid, 'ub',                  'num', pos(i), '- upper_bound')
    writeField(model, fid, 'grRules',             'txt', pos(i), '- gene_reaction_rule')
    writeField(model, fid, 'subSystems',          'txt', pos(i), '- subsystem')
    writeField(model, fid, 'rxnMiriams',          'txt', pos(i), '- annotation')
    writeField(model, fid, 'rxnConfidenceScores', 'num', pos(i), '- confidence_score')    
end

%Genes:
fprintf(fid,'- genes:\n');
[~,pos] = sort(model.genes);
for i = 1:length(model.genes)
    fprintf(fid,'  - !!omap\n');
    writeField(model, fid, 'genes',          'txt', pos(i), '- id')
    writeField(model, fid, 'geneShortNames', 'txt', pos(i), '- name')
end

%Compartments:
fprintf(fid,'- compartments:\n');
[~,pos] = sort(model.comps);
for i = 1:length(model.comps)
    writeField(model, fid, 'compNames', 'txt', pos(i), ['- ' model.comps{pos(i)}])
end

%TODO: include id, name & version (lost in RAVEN)

%Close file:
fclose(fid);

end

function writeField(model,fid,fieldName,type,pos,name)
%Writes a new line in the yaml file if the field exists and the field is
%not empty at the correspoinding position. It's recursive for some fields
%(metMiriams, rxnMiriams, and S)

if isfield(model,fieldName)
    if strcmp(fieldName,'metComps')
        %metComps: write full name
        fieldName = 'comps';
        pos       = model.metComps(pos);
    end
    
    field = eval(['model.' fieldName]);
    
    if strcmp(fieldName,'metMiriams')
        %metMiriams: create header & write chebi & kegg.compound
        if ~isempty(model.metCHEBI{pos}) || ~isempty(model.metKEGG{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            writeField(model, fid, 'metCHEBI', 'txt', pos, '  - chebi')
            writeField(model, fid, 'metKEGG',  'txt', pos, '  - kegg.compound')
        end
        
    elseif strcmp(fieldName,'rxnMiriams')
        %rxnMiriams: create header & write ec-codes & kegg.reaction
        if ~isempty(model.eccodes{pos}) || ~isempty(model.rxnKEGG{pos}) || ~isempty(model.rxnNotes{pos})
            fprintf(fid,['    ' name ': !!omap\n']);
            writeField(model, fid, 'eccodes',  'txt', pos, '  - ec-code')
            writeField(model, fid, 'rxnKEGG',  'txt', pos, '  - kegg.reaction')
            writeField(model, fid, 'rxnNotes', 'txt', pos, '  - pmid')
        end
        
    elseif strcmp(fieldName,'S')
        %S: create header & write each metabolite in a new line
        fprintf(fid,['    ' name ': !!omap\n']);
        if sum(field(:,pos) ~= 0) > 0
            model.mets   = model.mets(field(:,pos) ~= 0);
            model.coeffs = field(field(:,pos) ~= 0,pos);
            %Sort metabolites:
            [model.mets,order] = sort(model.mets);
            model.coeffs       = model.coeffs(order);
            for i = 1:length(model.mets)
                writeField(model, fid, 'coeffs',  'num', i, ['  - ' model.mets{i}])
            end
        end
        
    elseif sum(strcmp({'eccodes','rxnNotes','subSystems','metCHEBI','metKEGG','rxnKEGG'},fieldName)) > 0
        %eccodes/rxnNotes/subSystems: if 1 write in 1 line, if more create header and list
        if strcmp(fieldName,'subSystems')
            list = field{pos};  %subSystems already comes in a cell array
        else
            list = strrep(field{pos},' ','');     %Exception for eccodes & rxnNotes
            list = strsplit(list,';');
        end
        if length(list) == 1 && ~strcmp(list{1},'')
            fprintf(fid,['    ' name ': ' list{1} '\n']);
        elseif length(list) > 1
            fprintf(fid,['    ' name ':\n']);
            for i = 1:length(list)
                fprintf(fid,['        - ' list{i} '\n']);
            end
        end
        
    elseif sum(pos) > 0
        %All other fields:
        if strcmp(type,'txt')
            value = field{pos};
        elseif strcmp(type,'num')
            if isnan(field(pos))
                value = [];
            else
                value = num2str(field(pos));
            end
        end
        if ~isempty(value)
            fprintf(fid,['    ' name ': ' value '\n']);
        end
    end
end

end