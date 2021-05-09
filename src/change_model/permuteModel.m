function newModel=permuteModel(model, indexes, type)
% permuteModel
%   Changes the order of the reactions or metabolites in a model
%
%   Input:
%   model     a model structure
%   indexes   a vector with the same length as the number of items in the
%             model, which gives the new order of items
%   type      'rxns' for reactions, 'mets' for metabolites, 'genes' for
%             genes, 'comps' for compartments
%
% 	Output:
%   newModel  an updated model structure
%
% 	Usage: newModel=permuteModel(model, indexes, type)

newModel=model;
indexes=indexes(:);

switch type
    case 'rxns'
        if isfield(newModel,'rxns')
            newModel.rxns=newModel.rxns(indexes);
        end
        if isfield(newModel,'lb')
            newModel.lb=newModel.lb(indexes);
        end
        if isfield(newModel,'ub')
            newModel.ub=newModel.ub(indexes);
        end
        if isfield(newModel,'rev')
            newModel.rev=newModel.rev(indexes);
        end
        if isfield(newModel,'c')
            newModel.c=newModel.c(indexes);
        end
        if isfield(newModel,'S')
            newModel.S=newModel.S(:,indexes);
        end
        if isfield(newModel,'rxnNames')
            newModel.rxnNames=newModel.rxnNames(indexes);
        end
        if isfield(newModel,'rxnGeneMat')
            newModel.rxnGeneMat=newModel.rxnGeneMat(indexes,:);
        end
        if isfield(newModel,'grRules')
            newModel.grRules=newModel.grRules(indexes);
        end
        if isfield(newModel,'subSystems')
            newModel.subSystems=newModel.subSystems(indexes);
        end
        if isfield(newModel,'eccodes')
            newModel.eccodes=newModel.eccodes(indexes);
        end
        if isfield(newModel,'equations')
            newModel.equations=newModel.equations(indexes);
        end
        if isfield(newModel,'rxnMiriams')
            newModel.rxnMiriams=newModel.rxnMiriams(indexes);
        end
        if isfield(newModel,'rxnComps')
            newModel.rxnComps=newModel.rxnComps(indexes);
        end
        if isfield(newModel,'rxnFrom')
            newModel.rxnFrom=newModel.rxnFrom(indexes);
        end
        if isfield(newModel,'rxnScores')
            newModel.rxnScores=newModel.rxnScores(indexes);
        end
        if isfield(newModel,'rxnNotes')
            newModel.rxnNotes=newModel.rxnNotes(indexes);
        end
        if isfield(newModel,'rxnReferences')
            newModel.rxnReferences=newModel.rxnReferences(indexes);
        end
        if isfield(newModel,'rxnConfidenceScores')
            newModel.rxnConfidenceScores=newModel.rxnConfidenceScores(indexes);
        end
    case 'mets'
        if isfield(newModel,'mets')
            newModel.mets=newModel.mets(indexes);
        end
        if isfield(newModel,'metNames')
            newModel.metNames=newModel.metNames(indexes);
        end
        if isfield(newModel,'b')
            newModel.b=newModel.b(indexes,:);
        end
        if isfield(newModel,'metComps')
            newModel.metComps=newModel.metComps(indexes);
        end
        if isfield(newModel,'S')
            newModel.S=newModel.S(indexes,:);
        end
        if isfield(newModel,'unconstrained')
            newModel.unconstrained=newModel.unconstrained(indexes);
        end
        if isfield(newModel,'metMiriams')
            newModel.metMiriams=newModel.metMiriams(indexes,:);
        end
        if isfield(newModel,'inchis')
            newModel.inchis=newModel.inchis(indexes);
        end
        if isfield(newModel,'metFormulas')
            newModel.metFormulas=newModel.metFormulas(indexes);
        end
        if isfield(newModel,'metFrom')
            newModel.metFrom=newModel.metFrom(indexes);
        end
        if isfield(newModel,'metCharges')
            newModel.metCharges=newModel.metCharges(indexes);
        end
    case 'genes'
        if isfield(newModel,'genes')
            newModel.genes=newModel.genes(indexes);
        end
        if isfield(newModel,'geneComps')
            newModel.geneComps=newModel.geneComps(indexes);
        end
        if isfield(newModel,'geneMiriams')
            newModel.geneMiriams=newModel.geneMiriams(indexes);
        end
        if isfield(newModel,'geneShortNames')
            newModel.geneShortNames=newModel.geneShortNames(indexes);
        end
        if isfield(newModel,'rxnGeneMat')
            newModel.rxnGeneMat=newModel.rxnGeneMat(:,indexes);
        end
    case 'comps'
        if isfield(newModel,'comps')
            newModel.comps=newModel.comps(indexes);
        end
        if isfield(newModel,'compNames')
            newModel.compNames=newModel.compNames(indexes);
        end
        if isfield(newModel,'compOutside')
            newModel.compOutside=newModel.compOutside(indexes);
        end
        if isfield(newModel,'compMiriams')
            newModel.compMiriams=newModel.compMiriams(indexes);
        end
        [~,J]=sort(indexes); % The *index* of compartment is used in next fields
        if isfield(newModel,'metComps')
            [toreplace, bywhat] = ismember(newModel.metComps,1:length(J));
            newModel.metComps(toreplace) = J(bywhat(toreplace));
        end
        if isfield(model,'rxnComps')
            [toreplace, bywhat] = ismember(model.rxnComps,1:length(J));
            model.rxnComps(toreplace) = J(bywhat(toreplace));
        end
        if isfield(model,'geneComps')
            [toreplace, bywhat] = ismember(model.geneComps,1:length(J));
            model.geneComps(toreplace) = J(bywhat(toreplace));
        end
end
end
