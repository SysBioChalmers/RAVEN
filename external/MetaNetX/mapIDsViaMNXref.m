function subjectList=mapIDsViaMNXref(type,queryList,fromDB,toDB)
% mapIDsViaMNXref
%
%   Associate reaction/metabolite identifiers between different databases
%
%   type            'rxns' or 'mets' depending on which information is
%                   expected to associate
%   queryList       cell array of reactions/metabolites in the query database
%                   (e.g. model.rxnKEGGID, model.metChEBIID)
%   fromDB          the query database name
%   toDB            the subject database name
%
%   subjectList     cell array of reactions/metabolites in the subjuct database
%
%   Note: This function maps metabolite or reaction identifiers across
%   different databases using the MNXref naming convention. Each element in the
%   query cell array could be empty or contain only one identifier, which may
%   be mapped to multiple identifiers that are delimited by semicolon in the
%   subject cell array.
%
%   Usage: subjectList=mapIDsViaMNXref(type,queryList,fromDB,toDB)
%
%   Hao Wang, 2018-06-09
%

subjectList={};

if nargin<4
		EM='Missing input arguments';
		disp(EM);
end

if isequal(fromDB,toDB)
		EM='Query and subject databases cannot be the same!';
		disp(EM);
end

% Load MNXref data structure
load('MNXref.mat');

% Associate reaction or metabolite
if strcmpi(type,'rxns')
    MNXref=MNXrefRxns;
elseif strcmpi(type,'mets')
    MNXref=MNXrefMets;
else
    EM='Incorrect value of the "type" parameter. Allowed values are "rxns" or "mets"';
    dispEM(EM);
end

% Supported database names
dbNames=regexp(fieldnames(MNXref),'(\w+)MNXid','tokens');
dbNames=cellfun(@string,dbNames,'un',0);
dbNames=cellfun(@char,dbNames,'un',0);
dbNames=[dbNames;'MetaNetX'];
dbNames=setdiff(dbNames,{''});

if ~all(ismember({fromDB;toDB},dbNames))
		fprintf('Unknown database names! Please select from the following supported ones:\n');
		for i=1:numel(dbNames)
				fprintf('     %s\n',dbNames{i});
		end
		return;
end

%Prepare field names
fromDBField=strcat(fromDB,'xref');
fromMNXField=strcat(fromDB,'MNXid');
toMNXField=strcat(toDB,'MNXid');
toDBField=strcat(toDB,'xref');

% Initilize output cell array
subjectList=cell(numel(queryList),1);
subjectList(:)={''};

% In case of using MNX ids as query
if isequal('MetaNetX',fromDB)
		[a, b]=ismember(queryList,MNXref.(toMNXField));
		subjectList(find(a))=MNXref.(toDBField)(b(find(a)));
		dispNum(subjectList,MNXref.Version);
		return;
end

% Get the intermedia MNX identifiers
MNXids=cell(numel(queryList),1);
MNXids(:)={''};
[a, b]=ismember(queryList,MNXref.(fromDBField));
MNXids(find(a))=MNXref.(fromMNXField)(b(find(a)));

% In case of targting for MNX ids
if isequal('MetaNetX',toDB)
		subjectList=MNXids;
		dispNum(subjectList,MNXref.Version);
		return;
end

% Map intermedia MNX identifiers one by one
for i=1:numel(MNXids)
		if ~isempty(MNXids{i})
				index=strcmp(MNXids{i},MNXref.(toMNXField));
				if length(find(index))==1
						subjectList{i}=MNXref.(toDBField){find(index)};
				elseif length(find(index))>1
						subjectList{i}=strjoin(MNXref.(toDBField)(find(index)),';');
				end
		end
end
dispNum(subjectList,MNXref.Version);

end

% This subfunction counts and prints out the associated number
function dispNum(List,ver)
		num=numel(find(~cellfun(@isempty,List)));
		fprintf('%s ids were associated by MetaNetX version %s.\n',num2str(num),ver);
end
