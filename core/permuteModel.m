function newModel=permuteModel(model, indexes, type)
% permuteModel
%   Changes the order of the reactions or metabolites in a model
%
%   model     a model structure
%   indexes   a vector with the same length as the number of reactions in the
%             model which gives the new order of reactions
%   type      'rxns' for reactions and 'mets' for metabolites
%
% 	newModel  an updated model structure
%
% 	Usage: newModel=permuteModel(model, indexes, type)
%
%   Rasmus Agren, 2013-08-01
%

newModel=model;
indexes=indexes(:);

if strcmp(type,'rxns')
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
end

if strcmp(type,'mets')
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
end
end
