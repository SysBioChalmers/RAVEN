function checkFunctionUniqueness()
% checkFunctionUniqueness
%   Checks whether RAVEN toolbox functions are unique across other
%   toolboxes or user-created functions accessible in Matlab pathlist
%
%   Usage: checkFunctionUniqueness()

%Get the RAVEN path
ravenDir=findRAVENroot();

%Now getting all RAVEN functions recursively;
temp_res1=dir([ravenDir '/*/*.m']);
temp_res2=dir([ravenDir '/*/*/*.m']);

ravenFunctions={temp_res1.name,temp_res2.name}';
%startup.m is not a normal function, any startup.m in the path should run
%during startup, so duplicate use of this name is fine
ravenFunctions=ravenFunctions(~ismember(ravenFunctions,'startup.m'));

%Getting all the paths added to Matlab
if ispc
    matlabPaths=regexp(path, ';', 'split')';
elseif isunix
    matlabPaths=regexp(path, ':', 'split')';
end

hasConflicts=false;

for i=1:numel(matlabPaths)
    if ~any(strfind(matlabPaths{i},ravenDir))
        temp_res=dir([matlabPaths{i} '/*.m']);
        if ~isempty(temp_res)
            pathFunctions={temp_res.name}';
        else
            pathFunctions='';
        end
        if ~isempty(pathFunctions) && ~any(ismember('Contents.m',pathFunctions))
            if any(ismember(ravenFunctions,pathFunctions))
                fprintf('Not OK\n');
                if sum(ismember(ravenFunctions,pathFunctions))>(numel(ravenFunctions)/4)
                    EM='Multiple RAVEN versions detected in MATLAB path. Leave only one RAVEN version in MATLAB path and re-run checkInstallation\n';
                    dispEM(EM);
                else
                    disp(['WARNING: Duplicate functions in ',matlabPaths{i},': ']);
                    ovrlpFunctions=ravenFunctions(ismember(ravenFunctions,pathFunctions));
                    disp(ovrlpFunctions);
                    hasConflicts=true;
                end
            end
        end
    end
end

if hasConflicts
    fprintf('Pass\n')
    fprintf('  It is strongly recommended to resolve conflicting functions as this may compromise RAVEN functionality\n');
else
    fprintf('Fail\n');
end

end
