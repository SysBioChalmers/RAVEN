function steps = getINITSteps(metsToIgnore, series)
% getINITSteps
%   Converts a reaction score to the gene expression (CPM or TPM) required 
%   to get that reaction score, if the GPR is only a single gene.
%   Useful function primarily in test cases, where you want to be able to
%   define the reaction scores of rxns, but need to send in gene expression.
%   Note that all combinations of steps will not work. In general, avoid 'exclude'
%   if you want to define new ways to run the algorithm.
%
%   metsToIgnore  Structure describing mets that can be removed from the model
%                 before running ftINIT, such as water etc.
%                 (optional, default [])
%       simpleMets
%           mets  Names of metabolites to remove
%           compsToKeep Compartments for which metabolites should be kept.
%   series        Describes the way to run ftINIT: 
%                 '1+1'          Standard behavior. Step 1 and 2 described in
%                                the paper are merged into 1.
%                 '2+1'          The 3-step procedure described in the paper.
%                                Faster and slightly less accurate than '1+1 steps'
%                 '1+0'          Same as '1+1 steps', but skips step 3 in described 
%                                in the paper. This will result in a model including 
%                                a lot of reactions without GPRs. It is particularly 
%                                useful for structural comparison, since the reactions 
%                                removed in step 3 may be a bit random and doesn't 
%                                really add any extra information. Faster than
%                                '1+1 steps'
%                 '2+0'          Same as '2+1 steps', but skips step 3 in described 
%                                in the paper. Faster and slightly less accurate 
%                                than '1+0 steps', but will yield similar results.
%                 'full'         1-step run - similar to the old tINIT version but 
%                                without simplifications. Accurate, but very slow.
%                                This is mainly used for testing purposes.
%                 (optional, default '1+1')
%
%   steps         Cell array of steps, used as input to ftINIT
%
% Usage: steps = getINITSteps(metsToIgnore, series)
if nargin < 1
    metsToIgnore = [];
end

if nargin < 2
    series = '1+1';
end

if strcmp(series,'1+1') %step 1 and 2 are joined
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 120;
    params2 = struct();
    params2.MIPGap = 0.0030;
    params2.TimeLimit = 5000;
    params = {params1;params2};
    params3 = struct();
    params3.MIPGap = 0.0004;
    params3.TimeLimit = 5;
    paramsStep3 = {params3;params1;params2};
    %The paramsStep3 involves a quick first run. The objective value is often
    %small in the third step (~800), and 0.0004 of that is a very small number
    %With this first step, the rough value of the objective function will be 
    %estimated, which will generate an absolute MIPGap limit that is much larger
    %for the second iteration.
    
    steps = { ...
        INITStepDesc(false, false, 'ignore', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'essential', [1,0,0,0,1,0,0,0], metsToIgnore, paramsStep3, {10;10;20}); ...
        };
elseif strcmp(series,'2+1') %the 3 step process described in the ftINIT paper
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 120;
    params2 = struct();
    params2.MIPGap = 0.0030;
    params2.TimeLimit = 5000;
    params = {params1;params2};
    params3 = struct();
    params3.MIPGap = 0.0004;
    params3.TimeLimit = 5;
    paramsStep3 = {params3;params1;params2};
    %The paramsStep3 involves a quick first run. The objective value is often
    %small in the third step (~800), and 0.0004 of that is a very small number
    %With this first step, the rough value of the objective function will be 
    %estimated, which will generate an absolute MIPGap limit that is much larger
    %for the second iteration.
    
    steps = { ...
        INITStepDesc(true, true, 'ignore', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'essential', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'essential', [1,0,0,0,1,0,0,0], metsToIgnore, paramsStep3, {10;10;20}); ...
        };
elseif strcmp(series,'1+0') %Joins step 1 and 2, skips step 3
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 120;
    params2 = struct();
    params2.MIPGap = 0.0030;
    params2.TimeLimit = 5000;
    params = {params1;params2};
    steps = { ...
        INITStepDesc(false, false, 'ignore', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        };
elseif strcmp(series,'2+0') %Skips step 3
    params1 = struct();
    params1.MIPGap = 0.0004;
    params1.TimeLimit = 120;
    params2 = struct();
    params2.MIPGap = 0.0030;
    params2.TimeLimit = 5000;
    params = {params1;params2};
    steps = { ...
        INITStepDesc(true, true, 'ignore', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
        INITStepDesc(false, false, 'essential', [1,1,1,1,1,1,1,0], metsToIgnore, params, {10;20}); ...
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
