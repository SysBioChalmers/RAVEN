function modelSBML = constructYaml(modelSBML)
% constructYaml
%   Reduces the size of modelSBML structure before exporting to YAML format.
%   Removes empty and unnecessary fields.
%
%   modelSBML           model structure as generated during exportModel function
%
%   Usage: modelSBML = constructYaml(modelSBML)
%
%   Eduard Kerkhoven, 2018-03-19
%
useLessFields = {'spatialDimensions', 'size', 'isSetSize', ...
    'isSetSpatialDimensions', 'cvterms', 'initialAmount', ...
    'initialConcentration', 'substanceUnits', 'hasOnlySubstanceUnits', ...
    'constant', 'sboTerm', 'isSetInitialAmount', ...
    'isSetInitialConcentration', 'conversionFactor', 'level', ...
    'version', 'modifier', 'kineticLaw', 'fast', 'isSetFast', 'math', ...
    'constant', 'isSetValue', 'message'};

fnm = fieldnames(modelSBML);
idx = structfun(@isempty, modelSBML);
modelSBML = rmfield(modelSBML, fnm(idx));
fnm = fieldnames(modelSBML);
ws = modelSBML; % Keep working structure, as removing fields below will interfer with indexing

for i = 1:length(fnm)
    if isstruct(ws.(fnm{i}))
        fnm2 = fieldnames(ws.(fnm{i}));
        for j = 1:length(fnm2)
            if any(strcmp(fnm2{j}, useLessFields)) || ...
                    isempty([ws.(fnm{i})(:).(fnm2{j})])
                modelSBML.(fnm{i})=rmfield(modelSBML.(fnm{i}),fnm2{j});
%             elseif strcmp(fnm2{j}, 'fbc_geneProductAssociation')
%                 for l=1:length({ws.(fnm{i}).(fnm2{j})})
%                     modelSBML.(fnm{i})(l).(fnm2{j}).fbc_association = ...  
            elseif isstruct([ws.(fnm{i})(:).(fnm2{j})])
                for l=1:length({ws.(fnm{i}).(fnm2{j})})
                    if isstruct(ws.(fnm{i})(l).(fnm2{j}))
                        fnm3=fieldnames(ws.(fnm{i})(l).(fnm2{j}));
                        for k=1:length(fnm3) 
                            if any(strcmp(fnm3{k}, useLessFields)) || ...
                                    isempty([ws.(fnm{i})(l).(fnm2{j})(:).(fnm3{k})])
                                modelSBML.(fnm{i})(l).(fnm2{j}) = ...
                                    rmfield(modelSBML.(fnm{i})(l).(fnm2{j}),fnm3{k});
                            end
                        end
                    end
                end
            end
        end
    end
end
end