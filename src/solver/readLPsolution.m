function res = readLPsolution(solFile,logFile,solver)
% readLPsolution
%   Reads the SOL and LOG files written by the solver.

if nargin<3
    solver = 'gurobi';
end
if nargin<2
    logFile = '_tmp.log';
end
if nargin<1
    solFile = '_tmp.sol';
end

%Read .log file, different format for each solver
try
    fID = fopen(logFile);
    logData = textscan(fID,'%s','Delimiter', '\t');
    fclose(fID);
    logData = logData{1};
catch
    error(['Unable to load solver log file ' logFile])
end

switch solver
    case 'gurobi'
        %TODO: parse error messages
        obj = regexprep(logData{end-1}, 'Optimal objective ','');
        obj = str2double(obj);
end

%Read .sol file
try
    fID = fopen(solFile);
    solData = textscan(fID,'%s %f', 'Delimiter', ' ', 'CommentStyle', '#');
    fclose(fID);
    rxns = solData{1,1};
    flux = solData{1,2};
catch
    error(['Unable to load solver solution file ' solFile])
end

%x_.... are the reactions
keepIdx = startsWithOct(rxns,'x_');
rxns    = rxns(keepIdx);
flux    = flux(keepIdx);
rxns    = regexprep(rxns,'x_','');
rxns    = str2double(rxns);
[~,I]   = sort(rxns);
rxns    = rxns(I);
outFlux(rxns) = flux(I);
outFlux(~rxns) = 0;

res.full = outFlux;
res.obj  = obj;
delete(logFile);
delete(solFile);
