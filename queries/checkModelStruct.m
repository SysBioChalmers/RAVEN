function issues=checkModelStruct(model,varargin)
% checkModelStruct  Perform a number of checks to ensure a model structure is ok.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Name-Value Arguments
% --------------------
% throwErrors : logical
%     true if structural problems (wrong types, missing required fields,
%     duplicate IDs, etc.) should throw an error. Advisory notices such as
%     unused elements and bound inconsistencies are always reported as
%     warnings regardless of this flag (default true).
% trimWarnings : logical
%     true if only a maximum of 10 items should be displayed in a given
%     warning or error (default true).
%
% Returns
% -------
% issues : struct
%     When called with an output argument, returns a struct array with one
%     element per finding (does not throw or print) instead of the default
%     print/throw behaviour. Fields:
%
%     - category : char — type of issue: "missing_field", "wrong_type",
%       "empty_id", "duplicate", "invalid_id", "invalid_bounds",
%       "unused", "objective", "gpr", "invalid_formula",
%       "cross_reference", or "other".
%     - target : char — the specific item involved (field name, reaction
%       ID, metabolite ID, etc.), or "" when not applicable.
%     - message : char — the full diagnostic text.
%
% Notes
% -----
% This is typically performed after importing or constructing a model, or
% before attempting to export a model to SBML format.
%
% Examples
% --------
%     checkModelStruct(model);
%     checkModelStruct(model, "throwErrors", false);
%     issues = checkModelStruct(model);

p=parseRAVENargs(varargin, {'throwErrors',true; 'trimWarnings',true});
throwErrors=p.throwErrors;
trimWarnings=p.trimWarnings;

collecting = nargout > 0;
issues = struct('category',{},'target',{},'message',{});

    function reportIssue(severity, msg, items)
        % Report one issue: accumulate in collect mode, or warn/error.
        % severity: "error" | "warning"
        % When items is supplied (3 args) but is empty, silently do nothing.
        listProvided = nargin >= 3;
        if nargin < 3; items = {}; end
        if ~isempty(items)
            items = convertCharArray(items);
        end
        if listProvided && isempty(items); return; end
        if trimWarnings && numel(items) > 10
            items{10} = sprintf('...and %d more', numel(items)-9);
            items(11:end) = [];
        end
        cat = issueCategory(msg);
        if collecting
            if isempty(items)
                issues(end+1,1) = struct('category',cat,'target','','message',msg);
            else
                for ki_ = 1:numel(items)
                    issues(end+1,1) = struct('category',cat,'target',items{ki_},'message',msg); %#ok<AGROW>
                end
            end
            return
        end
        isErr = strcmp(severity,'error') && throwErrors;
        if isErr
            if isempty(items)
                error('RAVEN:modelError', '%s', msg);
            else
                error('RAVEN:modelError', '%s', ravenList(msg, items, false));
            end
        else
            if isempty(items)
                warning('RAVEN:modelError', '%s', msg);
            else
                warning('RAVEN:modelError', '%s', ravenList(msg, items, false));
            end
        end
    end

%Missing elements — checked inline so missing fields do not cause cascading
%errors in the rest of the body
fields={'id';'name';'rxns';'mets';'S';'lb';'ub';'rev';'c';'b';'comps';'metComps'};
for i=1:numel(fields)
    if ~isfield(model,fields{i})
        EM=['The model is missing the "' fields{i} '" field'];
        reportIssue('error',EM);
    end
end
if collecting && ~isempty(issues)
    return  % missing required fields would cause access errors below
end

%Type check
if ~ischar(model.id)
    EM='The "id" field must be a string';
    reportIssue('error',EM);
end
if ~ischar(model.name)
    EM='The "name" field must be a string';
    reportIssue('error',EM);
end
if ~iscellstr(model.rxns)
    EM='The "rxns" field must be a cell array of strings';
    reportIssue('error',EM);
end
if ~iscellstr(model.mets)
    EM='The "mets" field must be a cell array of strings';
    reportIssue('error',EM);
end
if ~isnumeric(model.S)
    EM='The "S" field must be of type "double"';
    reportIssue('error',EM);
end
if ~isnumeric(model.lb)
    EM='The "lb" field must be of type "double"';
    reportIssue('error',EM);
end
if ~isnumeric(model.ub)
    EM='The "ub" field must be of type "double"';
    reportIssue('error',EM);
end
if ~isnumeric(model.rev)
    EM='The "rev" field must be of type "double"';
    reportIssue('error',EM);
end
if ~isnumeric(model.c)
    EM='The "c" field must be of type "double"';
    reportIssue('error',EM);
end
if ~isnumeric(model.b)
    EM='The "b" field must be of type "double"';
    reportIssue('error',EM);
end
if ~iscellstr(model.comps)
    EM='The "comps" field must be a cell array of strings';
    reportIssue('error',EM);
end
if ~isnumeric(model.metComps)
    EM='The "metComps" field must be of type "double"';
    reportIssue('error',EM);
end
if isfield(model,'compNames')
    if ~iscellstr(model.compNames)
        EM='The "compNames" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'compOutside')
    if ~iscellstr(model.compOutside)
        EM='The "compOutside" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnNames')
    if ~iscellstr(model.rxnNames)
        EM='The "rxnNames" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'metNames')
    if ~iscellstr(model.metNames)
        EM='The "metNames" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'genes')
    if ~iscellstr(model.genes)
        EM='The "genes" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnGeneMat')
    if ~isnumeric(model.rxnGeneMat)
        EM='The "rxnGeneMat" field must be of type "double"';
        reportIssue('error',EM);
    end
end
if isfield(model,'grRules')
    if ~iscellstr(model.grRules)
        EM='The "grRules" field must be a cell array of strings';
        reportIssue('error',EM);
    end
    if ~isfield(model,'genes')
        EM='If "grRules" field exists, the model should also contain a "genes" field';
        reportIssue('error',EM);
    else
        %Erroneous grRules that start/end with OR/AND
        EM='The following reaction(s) have grRules that start or end with ''OR'' or ''AND'':';
        reportIssue('error',EM,model.rxns(startsWith(model.grRules,{'or ','and '}) | endsWith(model.grRules,{' or',' and'})));
        %grRules that are not in genes field
        geneList = getGenesFromGrRules(model.grRules);
        geneList = setdiff(unique(geneList),model.genes);
        if ~isempty(geneList)
            problemGrRules = model.rxns(contains(model.grRules,geneList));
            problemGrRules = strjoin(problemGrRules(:),'; ');
            EM=['The reaction(s) "' problemGrRules '" contain the following genes in its "grRules" field, but these are not in the "genes" field:'];
            reportIssue('error',EM,geneList);
        end
    end
end
if isfield(model,'rxnComps')
    if ~isnumeric(model.rxnComps)
        EM='The "rxnComps" field must be of type "double"';
        reportIssue('error',EM);
    end
end
if isfield(model,'inchis')
    if ~iscellstr(model.inchis)
        EM='The "inchis" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'metSmiles')
    if ~iscellstr(model.metSmiles)
        EM='The "metSmiles" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'metFormulas')
    if ~iscellstr(model.metFormulas)
        EM='The "metFormulas" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'metCharges')
    if ~isnumeric(model.metCharges)
        EM='The "metCharges" field must be a double';
        reportIssue('error',EM);
    end
end
if isfield(model,'metDeltaG')
    if ~isnumeric(model.metDeltaG)
        EM='The "metDeltaG" field must be a double';
        reportIssue('error',EM);
    end
end
if isfield(model,'subSystems')
    isNested  = any(cellfun(@(x) iscell(x), model.subSystems));
    isCellStr = any(cellfun(@(x) ischar(x), model.subSystems));
    if ~xor(isNested,isCellStr)
        EM='The "subSystems" field must be a cell array of chars, *or* a cell array of cell arrays of chars';
        reportIssue('error',EM);
    end
end
if isfield(model,'eccodes')
    if ~iscellstr(model.eccodes)
        EM='The "eccodes" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'unconstrained')
    if ~isnumeric(model.unconstrained)
        EM='The "unconstrained" field must be of type "double"';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnNotes')
    if ~iscellstr(model.rxnNotes)
        EM='The "rxnNotes" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnReferences')
    if ~iscellstr(model.rxnReferences)
        EM='The "rxnReferences" field must be a cell array of strings';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnConfidenceScores')
    if ~isnumeric(model.rxnConfidenceScores)
        EM='The "rxnConfidenceScores" field must be a double';
        reportIssue('error',EM);
    end
end
if isfield(model,'rxnDeltaG')
    if ~isnumeric(model.rxnDeltaG)
        EM='The "rxnDeltaG" field must be a double';
        reportIssue('error',EM);
    end
end

%Empty strings
if isempty(model.id)
    EM='The "id" field cannot be empty';
    reportIssue('error',EM);
end
if any(cellfun(@isempty,model.rxns))
    EM='The model contains empty reaction IDs';
    reportIssue('error',EM);
end
if any(cellfun(@isempty,model.mets))
    EM='The model contains empty metabolite IDs';
    reportIssue('error',EM);
end
if any(cellfun(@isempty,model.comps))
    EM='The model contains empty compartment IDs';
    reportIssue('error',EM);
end
if isfield(model,'metNames')
    EM='The following metabolites have empty names:';
    reportIssue('error',EM,model.mets(cellfun(@isempty,model.metNames)));
end

if isfield(model,'genes')
    if any(cellfun(@isempty,model.genes))
        EM='The model contains empty gene IDs';
        reportIssue('error',EM);
    end
end

%Validate format of ids
fields      = {'rxns';'mets';'comps';'genes'};
fieldNames  = {'reaction';'metabolite';'compartment';'gene'};
fieldPrefix = {'R_';'M_';'C_';'G_'};
for i=1:numel(fields)
    try
        numIDs = ~startsWith(model.(fields{i}),regexpPattern('^[a-zA-Z_]'));
    catch
        numIDs = [];
    end
    if any(numIDs)
        EM = ['The following ' fieldNames{i} ' identifiers do not start '...
            'with a letter or _ (conflicting with SBML specifications). '...
            'This does not impact RAVEN functionality, but be aware that '...
            'exportModel will automatically add ' fieldPrefix{i} ...
            ' prefixes to all ' fieldNames{i} ' identifiers:'];
        reportIssue('warning',EM,{model.(fields{i}){numIDs}});
    end
end

%Duplicates
EM='The following reaction IDs are duplicates:';
reportIssue('error',EM,model.rxns(duplicates(model.rxns)));
EM='The following metabolite IDs are duplicates:';
reportIssue('error',EM,model.mets(duplicates(model.mets)));
EM='The following compartment IDs are duplicates:';
reportIssue('error',EM,model.comps(duplicates(model.comps)));
if isfield(model,'genes')
    EM='The following genes are duplicates:';
    reportIssue('error',EM,model.genes(duplicates(model.genes)));
end
if isfield(model,'metNames')
    metInComp=strcat(model.metNames,'[',model.comps(model.metComps),']');
    EM='The following metabolites already exist in the same compartment:';
    reportIssue('error',EM,metInComp(duplicates(metInComp)));
end

%Elements never used (always reported as warnings)
EM='The following reactions are empty (no involved metabolites):';
reportIssue('warning',EM,model.rxns(~any(model.S,1)));
EM='The following metabolites are never used in a reaction:';
reportIssue('warning',EM,model.mets(~any(model.S,2)));
if isfield(model,'genes')
    EM='The following genes are not associated to a reaction:';
    reportIssue('warning',EM,model.genes(~any(model.rxnGeneMat,1)));
end
I=true(numel(model.comps),1);
I(model.metComps)=false;
EM='The following compartments contain no metabolites:';
reportIssue('warning',EM,model.comps(I));

%Contradicting bounds
EM='The following reactions have contradicting bounds (lower bound is higher than upper bound):';
reportIssue('error',EM,model.rxns(model.lb>model.ub));
EM='The following reactions have lower and upper bounds that indicate reversibility, but are indicated as irreversible in model.rev:';
reportIssue('warning',EM,model.rxns(model.lb < 0 & model.ub > 0 & model.rev==0));

%Multiple or no objective functions not allowed in SBML L3V1 FBCv2
if numel(find(model.c))>1
    EM='Multiple objective functions found. This might be intended, but results in FBCv2 non-compliant SBML file when exported';
    reportIssue('warning',EM,model.rxns(find(model.c)));
elseif ~any(model.c)
    EM='No objective function found. This might be intended, but results in FBCv2 non-compliant SBML file when exported';
    reportIssue('warning',EM);
end

%Mapping of compartments
if isfield(model,'compOutside')
    EM='The following compartments are in "compOutside" but not in "comps":';
    reportIssue('error',EM,setdiff(model.compOutside,[{''};model.comps]));
end

%Met names which start with number
if isfield(model,'metNames')
    I=false(numel(model.metNames),1);
    for i=1:numel(model.metNames)
        index=strfind(model.metNames{i},' ');
        if any(index)
            if any(str2double(model.metNames{i}(1:index(1)-1)))
                I(i)=true;
            end
        end
    end
    EM='The following metabolite names begin with a number directly followed by space, which could potentially cause problems:';
    reportIssue('warning',EM,model.metNames(I));
end

%Non-parseable composition
if isfield(model,'metFormulas')
    [~, ~, exitFlag]=parseFormulas(model.metFormulas,true,false);
    EM='The composition for the following metabolites could not be parsed:';
    reportIssue('warning',EM,model.mets(exitFlag==-1));
end

%Check if there are metabolites with different names but the same MIRIAM
%codes
if isfield(model,'metMiriams')
    miriams=containers.Map();
    for i=1:numel(model.mets)
        if ~isempty(model.metMiriams{i})
            %Loop through and add for each miriam
            for j=1:numel(model.metMiriams{i}.name)
                %Get existing metabolite indexes
                current=strcat(model.metMiriams{i}.name{j},'/',model.metMiriams{i}.value{j});
                if isKey(miriams,current)
                    existing=miriams(current);
                else
                    existing=[];
                end
                miriams(current)=[existing;i];
            end
        end
    end

    %Get all keys
    allMiriams=keys(miriams);

    hasMultiple=false(numel(allMiriams),1);
    for i=1:numel(allMiriams)
        if numel(miriams(allMiriams{i}))>1
            %Check if they all have the same name
            if numel(unique(model.metNames(miriams(allMiriams{i}))))>1
                if ~regexp(allMiriams{i},'^sbo\/SBO:') % SBO terms are expected to be multiple
                    hasMultiple(i)=true;
                end
            end
        end
    end

    %Print output
    EM='The following MIRIAM strings are associated to more than one unique metabolite name:';
    reportIssue('warning',EM,allMiriams(hasMultiple));
end

%Check if there are metabolites with different names but the same InChI
%codes
if isfield(model,'inchis')
    inchis=containers.Map();
    for i=1:numel(model.mets)
        if ~isempty(model.inchis{i})
            %Get existing metabolite indexes
            if isKey(inchis,model.inchis{i})
                existing=inchis(model.inchis{i});
            else
                existing=[];
            end
            inchis(model.inchis{i})=[existing;i];
        end
    end

    %Get all keys
    allInchis=keys(inchis);

    hasMultiple=false(numel(allInchis),1);
    for i=1:numel(allInchis)
        if numel(inchis(allInchis{i}))>1
            %Check if they all have the same name
            if numel(unique(model.metNames(inchis(allInchis{i}))))>1
                hasMultiple(i)=true;
            end
        end
    end

    %Print output
    EM='The following InChI strings are associated to more than one unique metabolite name:';
    reportIssue('warning',EM,allInchis(hasMultiple));
end

% %Check if there are metabolites with different names but the same SMILES
% if isfield(model,"metSmiles")
%     metSmiles=containers.Map();
%     for i=1:numel(model.mets)
%         if ~isempty(model.metSmiles{i})
%             %Get existing metabolite indexes
%             if isKey(metSmiles,model.metSmiles{i})
%                 existing=metSmiles(model.metSmiles{i});
%             else
%                 existing=[];
%             end
%             metSmiles(model.metSmiles{i})=[existing;i];
%         end
%     end
%
%     %Get all keys
%     allmetSmiles=keys(metSmiles);
%
%     hasMultiple=false(numel(metSmiles),1);
%     for i=1:numel(metSmiles)
%         if numel(metSmiles(metSmiles{i}))>1
%             %Check if they all have the same name
%             if numel(unique(model.metNames(metSmiles(allmetSmiles{i}))))>1
%                 hasMultiple(i)=true;
%             end
%         end
%     end
%
%     %Print output
%     EM="The following metSmiles strings are associated to more than one unique metabolite name:";
%     reportIssue("warning",EM,allmetSmiles(hasMultiple));
% end
end

function I=duplicates(strings)
I=false(numel(strings),1);
[J, K]=unique(strings);
if numel(J)~=numel(strings)
    L=1:numel(strings);
    L(K)=[];
    I(L)=true;
end
end

function cat=issueCategory(msg)
if contains(msg,'missing')
    cat='missing_field';
elseif contains(msg,'must be') || contains(msg,'should')
    cat='wrong_type';
elseif contains(msg,'empty')
    cat='empty_id';
elseif contains(msg,'duplicate') || contains(msg,'already exist in the same')
    cat='duplicate';
elseif contains(msg,'bound')
    cat='invalid_bounds';
elseif contains(msg,'SBML') || contains(msg,'identifier') || contains(msg,'do not start with')
    cat='invalid_id';
elseif contains(msg,'never used') || contains(msg,'not associated') || contains(msg,'contain no metabolites') || contains(msg,'are empty (no')
    cat='unused';
elseif contains(msg,'objective')
    cat='objective';
elseif contains(msg,'grRule') || contains(msg,'start or end with') || contains(msg,'"genes"')
    cat='gpr';
elseif contains(msg,'could not be parsed') || contains(msg,'composition')
    cat='invalid_formula';
elseif contains(msg,'MIRIAM') || contains(msg,'InChI') || contains(msg,'SMILES')
    cat='cross_reference';
else
    cat='other';
end
end
