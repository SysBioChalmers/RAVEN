function model=setParam(model, paramType, rxnList, params, var)
% setParam
%   Sets parameters for reactions
%
%   model       a model structure
%   paramType   the type of parameter to set:
%               'lb'    Lower bound
%               'ub'    Upper bound
%               'eq'    Both upper and lower bound (equality
%                       constraint)
%               'obj'   Objective coefficient
%               'rev'   Reversibility
%               'var'   Variance around measured bound
%   rxnList     a cell array of reaction IDs or a vector with their
%               corresponding indexes
%   params      a vector of the corresponding values
%   var         percentage of variance around measured value, if 'var' is
%               set as paramType. Defining 'var' as 5 results in lb and ub
%               at 97.5% and 102.5% of the provide params value (if params
%               value is negative, then lb and ub are 102.5% and 97.5%).
%
%   model       an updated model structure
%
%   Usage: model=setParam(model, paramType, rxnList, params)

if ~isempty(setdiff(paramType,{'lb';'ub';'eq';'obj';'rev';'var'}))
    EM=['Incorrect parameter type: "' paramType '"'];
    dispEM(EM);
end

%Allow to set several parameters to the same value
if numel(rxnList)~=numel(params) && numel(params)~=1
    EM='The number of parameter values and the number of reactions must be the same';
    dispEM(EM);
end

if isnumeric(rxnList) || islogical(rxnList)
    rxnList=model.rxns(rxnList);
elseif ischar(rxnList)
    rxnList={rxnList};
end

if ischar(paramType)
    paramType={paramType};
end

if isnumeric(params)
    params=[params];
end

if length(rxnList)>1
    if length(paramType)==1
        paramType(1:length(rxnList))=paramType;
    end
    if length(params)==1
        params(1:length(rxnList))=params;
    end
end

%Find the indexes for the reactions in rxnList. Do not use getIndexes
%as we do not want to throw errors if matches fail
indexes=zeros(numel(rxnList),1);

for i=1:numel(rxnList)
    index=find(strcmp(rxnList{i},model.rxns),1);
    if ~isempty(index)
        indexes(i)=index;
    else
        indexes(i)=-1;
        EM=['Reaction ' rxnList{i} ' is not present in the reaction list'];
        dispEM(EM,false);
    end
end

%Remove the reactions that were not found
params(indexes==-1)=[];
indexes(indexes==-1)=[];
paramType(indexes==-1)=[];
%Change the parameters

if ~isempty(indexes)
    if contains(paramType,'obj')
        model.c=zeros(numel(model.c),1); % parameter is changed, not added
    end     
    for j=1:length(paramType)
        if strcmpi(paramType{j},'eq')
            model.lb(indexes(j))=params(j);
            model.ub(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'lb')
            model.lb(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'ub')
            model.ub(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'obj')
            model.c(indexes(j))=params(j);
        end
        if strcmpi(paramType{j},'rev')
            %Non-zero values are interpreted as reversible
            model.rev(indexes(j))=params(j)~=0;
        end
        if strcmpi(paramType{j},'var')
            if params(j) < 0
                model.lb(indexes(j)) = params(j) * (1+var/200);
                model.ub(indexes(j)) = params(j) * (1-var/200);
            else
                model.lb(indexes(j)) = params(j) * (1-var/200);
                model.ub(indexes(j)) = params(j) * (1+var/200);
            end
        end
    end
end
end