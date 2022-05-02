function steps=getINITSteps(metsToIgnore, series)
if nargin < 1
    metsToIgnore = [];
end

if nargin < 2
    series = 'default';
end

if strcmp(series,'default') %the 3 step process described in the ftINIT paper
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 120;
    params2 = struct();
    params2.MIPGap = 0.0030;
    params2.TimeLimit = 5000;
    params = {params1;params2};
    steps = { ...
        INITStepDesc(true, true, 'ignore', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'exclude', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'essential', [1,0,0,0,1,0,0,0], metsToIgnore, params, {10;20}); ...
        };
elseif strcmp(series,'full') %Just one run, slow on large models, but this is the 'perfect' setup
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 10000;
    params = {params1};
    steps = { ...
        INITStepDesc(false, false, 'ignore', [0,0,0,0,0,0,0,0], [], params) ...
        };
    
else
    dispEM(['Invalid series in getINITSteps: ' series])
end

end
