function cModel=ravenToCobraModel(rModel)
% ravenToCobraModel
%   Converts RAVEN-compatible a constraint-based model structure to COBRA-compatible
%   structure
%
%   rModel        a RAVEN-compatible model structure
%
%   cModel        a COBRA-compatible model structure
%
%   The obtained COBRA model structure shouldn't be used on older COBRA
%   versions, which don't have FBCv2 compatibility. Such versions include
%   the commits older than April 15th 2016 in official COBRA for MATLAB
%   repository (https://github.com/opencobra/cobratoolbox/).
%
%   The fields are categorized like following (f.requiredEquiv):
%       I. Overlapping fields:
%           -rxns;
%           -mets;
%           -S;
%           -c;
%           -genes;
%           -lb;
%           -ub;
%           -rev;
%           -grRules;
%           -rxnGeneMat;
%           -metFormulas;
%           -description.
%       II. Fields which are not in RAVEN structure, but are mandatory in
%       COBRA structure (f.required):
%           -metCharge;
%           -rules;
%           -subSystems;
%       III. Optional COBRA fields with matching titles (f.optionalEquiv):
%           -id;
%           -rxnNames;
%           -metNames;
%           -rxnNotes;
%           -rxnReferences;
%           -confidenceScores.
%       IV. Optional fields with non-matching titles (f.optionalOtherName):
%           -metInchiString (inchis in RAVEN);
%           -rxnECNumbers (eccodes in RAVEN).
%
%   Usage: cModel = ravenToCobraModel(rModel)
%
%   Daniel Hermansson, 2016-08-29
%   Simonas Marcisauskas, 2017-01-31:
%       -fixed inconsistency with cModel.rules, when AND
%       relationships were ignored (function convertRules was completely
%       re-written);
%       -metabolites annotation information is now inherited in cModel;
%       -modified the whole function to make it compliant with additional
%       fields in RAVEN structure (accordingly to the newly added fields
%       rxnNotes, rxnReferences, confidenceScores and metCharge);
%       -removed the dependence from convFieldsCobra.m/confFieldsCobra.mat;
%       -made the function more comprehensive.
%

% Initializing f variable, which categorizes all the fields;
f.requiredEquiv={'rxns','mets','S','c','genes','lb','ub','rev','grRules','rxnGeneMat','metFormulas','description'};
f.required={'metCharge','rules','subSystems'};
f.optionalEquiv={'id','rxnNames','metNames','rxnNotes','rxnReferences','confidenceScores'};
f.optionalOtherName=containers.Map({'inchis','eccodes'},{'metInchiString','rxnECNumbers'});

% Creating a new COBRA-compatible structure, containing only the fields
% from group I, which are found in the current RAVEN model structure;
fieldOther = rmfield(rModel, intersect(fieldnames(rModel), f.requiredEquiv));
cModel = rmfield(rModel,fieldnames(fieldOther));

% Preparing cModel.mets..
% In COBRA structure model.mets must include the localisation information
% between square brackets, e.g. s_0001[c]. Concatenating metabolite ids and
% compartment abbreviations for that purpose;
cModel=setfield(cModel,'mets',convertMets(rModel.mets,rModel.comps(rModel.metComps)));

% Now preparing group II fields;
if (isfield(rModel,'metCharge')) 
	cModel=setfield(cModel,'metCharge',rModel.metCharge);
else
	cModel=setfield(cModel,'metCharge',zeros(numel(rModel.mets),1));
end
cModel=setfield(cModel,'rules',convertRules(rModel));
if (isfield(rModel,'subSystems')) 
	cModel=setfield(cModel,'subSystems',rModel.subSystems);
else
	cModel=setfield(cModel,'subSystems',repmat({''},size(rModel.rxns,1)));
end

% Now preparing group III fields;
fields = rmfield(rModel,setdiff(fieldnames(rModel), f.optionalEquiv));
cModel=structUpdate(cModel,fields);

% Now preparing group IV fields;
if (isfield(rModel,'inchis')) cModel=setfield(cModel,f.optionalOtherName('inchis'),rModel.inchis); end
	if (isfield(rModel,'eccodes')) cModel=setfield(cModel,f.optionalOtherName('eccodes'),rModel.eccodes); end

% Now obtaining the annotation for metabolites,
% which is stored in rModel.metMiriams field. From this field
% cModel.metCHEBIID and cModel.metKEGGID are be generated.
cModel.metKEGGID=getMiriamVector(rModel,'metKegg');
cModel.metCHEBIID=getMiriamVector(rModel,'metChebi');

% The final step is to get annotation for reaction references. Several
% pubmed referenced could be in rModel.rxnReferences already, so if there
% are any references in rModel.rxnMiriams, these must be concatenated with
% cModel.rxnReferences, which was obtain during group III fields
% processing;
if isfield(cModel,'rxnReferences')
    cModel.rxnReferences=strcat(cModel.rxnReferences,';',getMiriamVector(rModel,'rxnPubmed',false,'pubmed:'));
    cModel.rxnReferences=regexprep(cModel.rxnReferences,'^;','');
else
    cModel.rxnReferences=getMiriamVector(rModel,'rxnPubmed',false,'pubmed:');
end
end

function mets=convertMets(mets,metComps)
    % Removing the trailing abbreviation in mets about
    % compartmentalization, for instance, s_0001_c becomes s_0001. It is
    % assumed that abbreviation is no longer than 3 characters;
    mets=regexprep(mets,'_.{1,3}$','');
    % Now concatenating metabolite ids with compartmentalization
    % information, available in model.metComps;
	mets=arrayfun(@(x) [mets{x} '[' metComps{x} ']'],[1:numel(mets)],'UniformOutput',false);
	mets=transpose(mets);
end

function rules=convertRules(model)
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

function s_merged=structUpdate(s_old,s_new)
 	%// Remove overlapping fields from first struct%// Obtain all unique names of remaining fields,%// Merge both structs
 	s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
 	names = [fieldnames(s_merged); fieldnames(s_new)];
 	s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end