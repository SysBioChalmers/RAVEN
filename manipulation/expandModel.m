function [newModel, rxnToCheck]=expandModel(model)
% expandModel  Expand reactions that use several gene associations.
%
% Each reaction that uses several gene associations is split into several
% reactions, each under the control of only one gene.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Returns
% -------
% newModel : struct
%     model structure with separate reactions for iso-enzymes, where the
%     reaction ids are renamed as id_EXP_1, id_EXP_2, etc.
% rxnToCheck : cell
%     cell array with original reaction identifiers for those whose grRule
%     was not in disjunctive normal form, and so had to be expanded by
%     distributivity, plus any whose grRule could not be parsed.
%
% Examples
% --------
%     [newModel, rxnToCheck]=expandModel(model);
%
% Notes
% -----
% The isozymes are the AND-clauses of the grRule in disjunctive normal form,
% so nested expressions are expanded correctly: "g1 and (g2 or g3)" yields
% "g1 and g2" and "g1 and g3". Reactions listed in rxnToCheck are expanded
% correctly too; they are reported only because a rule needing
% distributivity is often a sign that the GPR was not what its author
% intended. Use standardizeGrRules/findPotentialErrors to inspect them.

%Work out the isozymes of every reaction up front. The number of copies is
%the number of DNF clauses, which is not the number of ' or ' substrings:
%"(g1 or g2) and (g3 or g4)" has two or:s but four isozymes.
prevNumRxns = length(model.rxns);
clauses = cell(prevNumRxns,1);
nCopies = zeros(prevNumRxns,1);
rxnToCheck={};
for i=1:prevNumRxns
    rule = model.grRules{i};
    if isempty(rule)
        continue
    end
    try
        c = grRuleToDNF(rule);
    catch ME
        if ~strcmp(ME.identifier,'RAVEN:badGrRule')
            rethrow(ME)
        end
        %An unparseable rule is left untouched rather than guessed at.
        rxnToCheck{end+1,1}=model.rxns{i}; %#ok<AGROW>
        continue
    end
    if numel(c) <= 1
        continue
    end
    %DEVIATION from raven-toolbox expand.py, which has no such guard: the
    %cross-product is exponential in the number of ORed complexes, and
    %MATLAB grows the model fields eagerly, so a pathological rule would
    %exhaust memory instead of raising. Fail with the culprit named.
    if numel(c) > 10000
        error('RAVEN:grRuleTooComplex',['Reaction ' model.rxns{i} ' has a grRule ' ...
            'with ' num2str(numel(c)) ' isozymes after expansion. This is almost ' ...
            'certainly a malformed rule; check it with findPotentialErrors.']);
    end
    clauses{i} = c;
    nCopies(i) = numel(c)-1;
    if ~isDnfGrRule(rule)
        rxnToCheck{end+1,1}=model.rxns{i}; %#ok<AGROW>
    end
end
toAdd = sum(nCopies);
if toAdd > 0
    %Calculate indices to expand
    %For example, if a reaction with index x has 2 or:s, meaning it has 3
    %reactions after the split, we should add two copies of this reaction
    %For fields that should just be copied to the new reactions, we just keep
    %track of that there are two copies, i.e., we add x x to this vector.
    %That is exactly what repelem does for us.
    %(:) forces a column. Indexing a 1x1 field takes the orientation of the
    %index rather than of the field, so a row here made every concatenation
    %below fail on a single-reaction model.
    cpyIndices = reshape(repelem((1:prevNumRxns)', nCopies), [], 1);
    
    %Copy all fields that should just be copied
    model.S=[model.S model.S(:,cpyIndices)];
    model.rxnNames=[model.rxnNames;model.rxnNames(cpyIndices)];
    model.lb=[model.lb;model.lb(cpyIndices)];
    model.ub=[model.ub;model.ub(cpyIndices)];
    model.rev=[model.rev;model.rev(cpyIndices)];
    model.c=[model.c;model.c(cpyIndices)];
    if isfield(model,'subSystems')
        model.subSystems=[model.subSystems;model.subSystems(cpyIndices)];
    end
    if isfield(model,'eccodes')
        model.eccodes=[model.eccodes;model.eccodes(cpyIndices)];
    end
    if isfield(model,'equations')
        model.equations=[model.equations;model.equations(cpyIndices)];
    end
    if isfield(model,'rxnMiriams')
        model.rxnMiriams=[model.rxnMiriams;model.rxnMiriams(cpyIndices)];
    end
    if isfield(model,'rxnComps')
        model.rxnComps=[model.rxnComps;model.rxnComps(cpyIndices)];
    end
    if isfield(model,'rxnFrom')
        model.rxnFrom=[model.rxnFrom;model.rxnFrom(cpyIndices)];
    end
    if isfield(model,'rxnNotes')
        model.rxnNotes=[model.rxnNotes;model.rxnNotes(cpyIndices)];
    end
    if isfield(model,'rxnReferences')
        model.rxnReferences=[model.rxnReferences;model.rxnReferences(cpyIndices)];
    end
    if isfield(model,'rxnConfidenceScores')
        model.rxnConfidenceScores=[model.rxnConfidenceScores;model.rxnConfidenceScores(cpyIndices)];
    end
    if isfield(model,'rxnDeltaG')
        model.rxnDeltaG=[model.rxnDeltaG;model.rxnDeltaG(cpyIndices)];
    end
    
    %now expand the more complex fields - will be filled in later
    model.rxns=[model.rxns;cell(toAdd,1)];
    model.grRules=[model.grRules;cell(toAdd,1)];
    model.rxnGeneMat=[model.rxnGeneMat;sparse(toAdd,size(model.rxnGeneMat,2))];
    
    %Loop throught those reactions and fill in the expanded data
    nextIndex = prevNumRxns + 1;
    for i=1:prevNumRxns
        if (nCopies(i) > 0)
            geneSets = clauses{i};

            %The first isozyme replaces the original reaction
            model.grRules{i}=strjoin(geneSets{1},' and ');
            model.rxnGeneMat(i,:)=0;
            model.rxnGeneMat(i,ismember(model.genes,geneSets{1})) = 1;

            %The rest are appended, named after the original reaction, which
            %is why it is only renamed once the loop is done
            for j=2:numel(geneSets)
                ind = nextIndex+j-2;
                model.rxns{ind}=[model.rxns{i} '_EXP_' num2str(j)];
                model.grRules{ind}=strjoin(geneSets{j},' and ');
                model.rxnGeneMat(ind,ismember(model.genes,geneSets{j})) = 1;
            end
            model.rxns{i}=[model.rxns{i}, '_EXP_1'];
            nextIndex = nextIndex + nCopies(i);
        end
    end
    newModel=model;
else
    %There are no reactions to expand, return the model as is
    newModel=model;
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(newModel,true);
newModel.grRules = grRules;
newModel.rxnGeneMat = rxnGeneMat;
end
