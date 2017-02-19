function metabolites=parseRxnEqu(equations)
% parseRxnEqu
%   Gets all metabolite names from a cell array of equations
%
%   metabolites=parseRxnEqu(equations)
%
%   equations     A cell array with equation strings
%
%   metabolites   A cell array with the involved metabolites
%
%   The equations should be written like:
%   1 A + 3 B (=> or <=>) 5C + 2 D
%
%   If the equation is expressed as for example '... + (n-1) starch' then
%   '(n-1) starch' will be interpreted as one metabolite
%
%   Usage: metabolites=parseRxnEqu(equations)
%
%   Rasmus Agren, 2012-05-27
%

if ~iscell(equations)
    equations={equations};
end

metabolites={};

%Replace the the direction arrows and plus signs with a weird character
%that will be used for parsing
equations=strrep(equations,' <=> ', '$$$');
equations=strrep(equations,' => ', '$$$');
equations=strrep(equations,' + ', '$$$');
equations=strtrim(equations);

for i=1:numel(equations)
    %Split each equation in possible metabolites
    candidates=regexp(equations{i},'$$$','split');

    %If the splitting character is at the end (if exchange rxns), then an
    %empty string will exist together with the real ones. Remove it
    candidates(cellfun(@isempty,candidates))=[];

    %Now remove the potential coefficient before each metabolite
    for j=1:numel(candidates)
        %If the metabolite has a coefficient it will look as 'number name'
        space=strfind(candidates{j},' ');

        if isempty(space)
            %Add the metabolite
            metabolites=[metabolites;candidates(j)];
        else
            potNumber=candidates{j}(1:space(1));
            %I use str2double here which can't deal with fractions (1/3
            %glc and so on). I do this because I don't want to risk calling
            %functions
            [~,isNumber]=str2num(potNumber);

            if isNumber==1
                %Remove the coefficient
                metName=candidates{j}(space(1)+1:end);
                metabolites=[metabolites;metName];
            else
                %The metabolite name contained spaces
                metabolites=[metabolites;candidates(j)];
            end
        end
    end
end

metabolites=strtrim(metabolites);

%Upper/lower case is treated as different names. This should be checked for
%later since it's bad modelling practice
metabolites=unique(metabolites);

end
