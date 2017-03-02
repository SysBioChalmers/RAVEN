function miriamVector=getMiriamVector(model,selection,addNull,header)
% getMiriamVector
%   Obtains a vector, which contains KEGG ids from rxnMiriams or
%   metMiriams. This function is useful if ones wants to unfold metMiriams
%   and rxnMiriams to the simple cell array, like metCHEBIID in COBRA
%   metabolic model structure
%
%   model             a model structure
%   selection         a string, which states, whether miriam vector should
%                     be provided for reactions or metabolites.
%                     Available options: 'rxnKegg', 'rxnPubmed, 'metKegg',
%                     'metChebi'
%   addNull           true if the entries for reactions / metabolites
%                     without any annotation data should be described as
%                     concatenated string of model.id and ':null', otherwise
%                     the corresponding entry is left blank (optional,
%                     default false)
%   header            the string, which precedes each miriam entry (optional, default '')
%
%   miriamVector      a vector with KEGG / CHEBI ids
%
%   Usage: miriamVector=getMiriamVector(model,selection,addNull,header)
%
%   Simonas Marcisauskas, 2015-09-23
%

if nargin<3
	addNull=false;
end;
if nargin<4
	header='';
end;

% Firstly creating rawHeaders, useful to locate the specific miriams;
if strcmp(selection,'rxnKegg')
    rawHeader='kegg.reaction';
elseif strcmp(selection,'rxnPubmed')
    rawHeader='pubmed';
elseif strcmp(selection,'metKegg')
    rawHeaderA='kegg.compound';
    rawHeaderB='kegg.glycan';
elseif strcmp(selection,'metChebi')
    rawHeaderA='obo.chebi:CHEBI';
    rawHeaderB='chebi:CHEBI';
else
    error('Miriam type is not defined in the function. See the function for more details');
end;

% The main part of the script;
if strcmp(selection,'rxnKegg') || strcmp(selection,'rxnPubmed')
    % Creating an empty cell array for rxn ids vector;
    % Firstly obtaining the list of relevant miriam ids. Several reactions
    % may have several miriam ids, such ids are kept in additional columns;
    rxnMiriams=cell([numel(model.rxns) 1]);
    for i=1:numel(model.rxns)
        if (~isempty(model.rxnMiriams{i,1})) && any(strcmp(model.rxnMiriams{i,1}.name,rawHeader))
            for j=1:numel(model.rxnMiriams{i,1}.name)
                if strcmp(model.rxnMiriams{i,1}.name(j),rawHeader)
                    rxnMiriams(i,j) = model.rxnMiriams{i,1}.value(j);
                end;           
            end;
        else
            rxnMiriams{i,1} = strcat(model.id,':null');
            if addNull==false
            	rxnMiriams{i,1} = regexprep(rxnMiriams{i,1},strcat(model.id,':null'),'');
            end;
        end;
    end;
    % Now adding headers to the newly obtained miriam ids. For reactions
    % which are associated with more than one miriam, these are
    % concatenated by using semicolon as separator;
    tempMiriams = cell([size(rxnMiriams,1) 1]);
    for i=1:size(rxnMiriams,1)
    	for j=1:size(rxnMiriams,2)
            if j==1
            	tempMiriams{i,1}=strcat(header,rxnMiriams{i,1},';');
            else
                tempMiriams{i,1}=strcat(tempMiriams{i,1},header,rxnMiriams{i,j},';');
            end;
        end;
        if nargin==4
            tempMiriams{i,1} = regexprep(tempMiriams{i,1},strcat(header,';|$'),'');
        end;
        tempMiriams{i,1} = regexprep(tempMiriams{i,1},'^;|;+$','');
    end;
    miriamVector=tempMiriams;
elseif strcmp(selection,'metKegg') || strcmp(selection,'metChebi')
% Creating an empty cell array for met ids vector;
% Firstly obtaining the list of relevant miriam ids. Several metabolites
% may have several miriam ids, such ids are kept in additional columns;
    metMiriams=cell([numel(model.mets) 1]);
    for i=1:numel(model.mets)
        if (~isempty(model.metMiriams{i,1})) && (any(strcmp(model.metMiriams{i,1}.name,rawHeaderA)) || any(strcmp(model.metMiriams{i,1}.name,rawHeaderB)))
            for j=1:numel(model.metMiriams{i,1}.name)
                if strcmp(model.metMiriams{i,1}.name(j),rawHeaderA) || strcmp(model.metMiriams{i,1}.name(j),rawHeaderB)
                    metMiriams(i,j) = model.metMiriams{i,1}.value(j);
                end;           
            end;
        else
            metMiriams{i,1} = strcat(model.id,':null');
            if addNull==false
                metMiriams{i,1} = regexprep(metMiriams{i,1},strcat(model.id,':null'),'');
            end;
        end;
    end;
    % Now adding headers to the newly obtained miriam ids. For metabolites
    % which are associated with more than one miriam, these are
    % concatenated by using semicolon as separator;
    tempMiriams = cell([size(metMiriams,1) 1]);
    for i=1:size(metMiriams,1)
        for j=1:size(metMiriams,2)
            if j==1
                tempMiriams{i,1}=strcat(header,metMiriams{i,1},';');
            else
                tempMiriams{i,1}=strcat(tempMiriams{i,1},header,metMiriams{i,j},';');
            end;
        end;
        if nargin==4
            tempMiriams{i,1} = regexprep(tempMiriams{i,1},strcat(header,';|$'),'');
        end;
        tempMiriams{i,1} = regexprep(tempMiriams{i,1},'^;|;+$','');
    end;
    miriamVector=tempMiriams;
end;
end