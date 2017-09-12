function checkFunctionUniqueness()
% checkFunctionUniqueness
%   Checks whether RAVEN toolbox functions are unique across other
%   toolboxes or user-created functions accessible in Matlab pathlist
%
%   Usage: checkFunctionUniqueness()
%
%   Simonas Marcisauskas, 2017-09-08
%

%Get the RAVEN path
[ST, I]=dbstack('-completenames');
[ravenDir,~,~]=fileparts(fileparts(ST(I).file));

%Now getting all RAVEN functions recursively;
temp_res1=dir([ravenDir '/*/*.m']);
temp_res2=dir([ravenDir '/*/*/*.m']);

ravenFunctions={temp_res1.name,temp_res2.name}';

%Getting all the paths added to Matlab
matlabPaths=regexp(path, ':', 'split')';

hasConflicts=false;

for i=1:numel(matlabPaths)
    if ~contains(matlabPaths{i},ravenDir)
        temp_res=dir([matlabPaths{i} '/*.m']);
        if ~isempty(temp_res)
            pathFunctions={temp_res.name}';
        else
            pathFunctions='';
        end;
        if ~isempty(pathFunctions)
            if any(ismember(ravenFunctions,pathFunctions))
                fprintf(['Duplicate functions in ',matlabPaths{i},':\n']);
                ovrlpFunctions=ravenFunctions(ismember(ravenFunctions,pathFunctions));
                disp(ovrlpFunctions);
                hasConflicts=true;
            end;
        end;
    end;
end;

if hasConflicts
    fprintf('It is strongly recommended to resolve conflicting functions as this may compromise RAVEN functionality\n');
end;

end
