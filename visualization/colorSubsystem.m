function newMap = colorSubsystem(map, model, subsystem, color)
% colorSubsystem  Color every reaction in a specific subsystem.
%
% Modified from colorSubsystemCD as distributed through the COBRA Toolbox
% (https://github.com/opencobra/cobratoolbox/blob/src/visualization/metabolicCartography/colorSubsystemCD.m).
%
% Parameters
% ----------
% map : struct
%     file from CellDesigner parsed to MATLAB format.
% model : struct
%     COBRA model structure.
% subsystem : char
%     name of a subsystem as a string.
% color : char, optional
%     color desired for reactions, in capitals (default 'RED').
%
% Returns
% -------
% newMap : struct
%     MATLAB structure of map with reaction modifications.
%
% Examples
% --------
%     newMap = colorSubsystem(map, model, subsystem, color);
if nargin < 4
    color = 'RED';
end

newMap = map;
rxnList = model.rxns(ismember([model.subSystems{:}]', subsystem));
colors = createColorsMap;

index = find(ismember(newMap.rxnName, rxnList));
for j = index'
    newMap.rxnColor{j, 1} = colors(color);
    %newMap.rxnWidth{j, 1} = areaWidth;
end

% Use the existence of reactant lines to check if the map has the
% complete structure, and if so change also secondary lines.
if any(strcmp('rxnReactantLineColor', fieldnames(newMap))) == 1
    for j = index'
        if ~isempty(newMap.rxnReactantLineColor{j})
            for k = 1:length(newMap.rxnReactantLineColor{j})
                newMap.rxnReactantLineColor{j, 1}{k, 1} = colors(color);
                %newMap.rxnReactantLineWidth{j, 1}{k, 1} = areaWidth;
            end
        end
        if ~isempty(newMap.rxnProductLineColor{j})
            for m = 1:1:length(newMap.rxnProductLineColor{j})
                newMap.rxnProductLineColor{j, 1}{m, 1} = colors(color);
                %newMap.rxnProductLineWidth{j, 1}{m, 1} = areaWidth;
            end
        end
    end
end

end
