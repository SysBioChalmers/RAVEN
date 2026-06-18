function newModel=permuteModel(model, indexes, type)
% permuteModel  Change the order of the reactions or metabolites in a model.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% indexes : double
%     a vector with the same length as the number of items in the model,
%     which gives the new order of items.
% type : char
%     'rxns' for reactions, 'mets' for metabolites, 'genes' for genes,
%     'comps' for compartments.
%
% Returns
% -------
% newModel : struct
%     an updated model structure.
%
% Examples
% --------
%     newModel = permuteModel(model, indexes, type);

newModel=model;
type=char(type);

reg=ravenModelFields();
regKey=type; if regKey(end)=='s', regKey=regKey(1:end-1); end
typeReg=reg(strcmp({reg.type},regKey));

% Apply the permutation to all 1-D registry fields of this type.
for i=1:numel(typeReg)
    fname=typeReg(i).name;
    if isfield(newModel,fname)
        newModel.(fname)=newModel.(fname)(indexes);
    end
end

% Handle 2-D and value-remapping fields not covered by the generic loop.
switch type
    case 'rxns'
        if isfield(newModel,'S')
            newModel.S=newModel.S(:,indexes);
        end
        if isfield(newModel,'rxnGeneMat')
            newModel.rxnGeneMat=newModel.rxnGeneMat(indexes,:);
        end
    case 'mets'
        if isfield(newModel,'S')
            newModel.S=newModel.S(indexes,:);
        end
        if isfield(newModel,'b')
            newModel.b=newModel.b(indexes,:);
        end
    case 'genes'
        if isfield(newModel,'rxnGeneMat')
            newModel.rxnGeneMat=newModel.rxnGeneMat(:,indexes);
        end
    case 'comps'
        % comps/compNames/compOutside/compMiriams were permuted above.
        % Remap VALUE-arrays that store indices into comps.
        [~,J]=sort(indexes);
        for fn={'metComps','rxnComps','geneComps'}
            fname=fn{1};
            if isfield(newModel,fname)
                [toreplace,bywhat]=ismember(newModel.(fname),1:length(J));
                newModel.(fname)(toreplace)=J(bywhat(toreplace));
            end
        end
end
end
