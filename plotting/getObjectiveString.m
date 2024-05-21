function objectiveString = getObjectiveString(max, objectiveNames, objectiveValues)
% getObjectiveString
%   Returns a string representing the objective function (e.g. 'MAX Growth
%   - 0.5 HXT4').
%
%   max                 true if the objective function should be maximized
%   objectiveNames      cell array of reaction names
%   objectiveValues     the corresponding coefficients for each reaction
%
%   objectiveString     the calculated objective function
%
% Usage: objectiveString = getObjectiveString(max, objectiveNames,
%           objectiveValues)

objectiveString='';

if max==true
    objectiveString=[objectiveString 'MAX: '];
else
    objectiveString=[objectiveString 'MIN: '];
end

%Loops through the reactions
for i=1:length(objectiveNames)
    %Add no sign if it's the first reaction
    if i>1
        if objectiveValues(i)==1
            objectiveString=[objectiveString ' + ' objectiveNames{i}];
            continue;
        end
        if objectiveValues(i)==-1
            objectiveString=[objectiveString ' - ' objectiveNames{i}];
            continue;
        end
        if objectiveValues(i)>=0
            objectiveString=[objectiveString ' + ' num2str(objectiveValues(i)) ' ' objectiveNames{i}];
        else
            objectiveString=[objectiveString ' - ' num2str(abs(objectiveValues(i))) ' ' objectiveNames{i}];
        end
    else
        if objectiveValues(i)==1
            objectiveString=[objectiveString objectiveNames{i}];
            continue;
        end
        if objectiveValues(i)==-1
            objectiveString=[objectiveString '- ' objectiveNames{i}];
            continue;
        end
        if objectiveValues(i)>=0
            objectiveString=[objectiveString num2str(objectiveValues(i)) ' ' objectiveNames{i}];
        else
            objectiveString=[objectiveString '- ' num2str(abs(objectiveValues(i))) ' ' objectiveNames{i}];
        end
    end
end
end
