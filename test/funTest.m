function out=funTest(funCell,parCell,logs)
% runTest
%  Run unit tests specified as functions.
%
%   funCell	cell containing the names of tests to be run
%   parCell	test specific parameters.
%   		Use an empty cell for no parameters.
%   logs	a boolean to specify if output should be logged
%			(default true). Logs are saved to './log/'.
%
%   outStruct	output struct containing test specific results
%
%   The test functions should be stored in folder 'private' in
%   the same direcory as this function.
%
%   Usage:
%   fun={'fun1' 'fun2'};
%   param={struct('someParam',13),'someParamForFun2'};
%   out=runTest(fun,param,true)
%
%   Daniel Hermansson, 2016-04-19

if nargin<3
    logs=true;
end
if nargin<2
    parCell=cell(1,length(funCell));
end

fileID=generateID(6); % generate a random log file identifier

out={};
pass={};
time=[];

for i=1:numel(parCell)
    if i<=length(funCell)
        fun=str2func(funCell{i});
        funstr=funCell{i};
    end
    
    if(logs)
        diary(['./log/',fileID,'_',funstr,'_',num2str(i),'.log']);
    end
    
    disp(['########## NEW TEST: ' funstr ' ##########'])
    fprintf(['Params:\n'])
    if (~isempty(parCell{i})) disp(parCell{i}); end
    
    try
        tic;
        out{i}=fun(parCell{i})
        time(i)=toc;
        pass{i}='PASS';
    catch EM
        time(i)=toc;
        fprintf([funstr,' failed execution. Params:\n'])
        disp(parCell{i})
        disp(EM)
        pass{i}='FAIL';
    end
    
    if(logs)
        diary off
        diary(['./log/',fileID,'_Summary.log']);
    end
    disp(['####################'])
    fprintf([funstr '_' num2str(i) ':	' pass{i} '\n'])
    fprintf(['Params:\n'])
    if (~isempty(parCell{i})) disp(parCell{i}); end
    fprintf(['Runtime:	' num2str(time(i)) ' s\n\n'])
end

diary off
end

% supporting functions
function R=generateID(N)
% mix numbers and characters -> use characters
SET = char(['a':'z' '0':'9']);
NSET = length(SET);

%pick N numbers
i = ceil(NSET*rand(1,N)); % with repeat
R = SET(i);
end

function s_merged=structUpdate(s_old,s_new)
s_merged = rmfield(s_old, intersect(fieldnames(s_old), fieldnames(s_new)));
names = [fieldnames(s_merged); fieldnames(s_new)];
s_merged = cell2struct([struct2cell(s_merged); struct2cell(s_new)], names, 1);
end
