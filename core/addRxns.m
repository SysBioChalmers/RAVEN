function newModel=addRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
% addRxns
%   Adds reactions to a model
%
%   model            a model structure
%   rxnsToAdd        the reaction structure can have the following fields:
%            rxns               cell array with unique strings that
%                               identifies each reaction
%            equations          cell array with equation strings. Decimal
%                               coefficients are expressed as "1.2".
%                               Reversibility is indicated by "<=>" or "=>"
%            rxnNames           cell array with the names of each reaction
%                               (opt, default '')
%            lb                 vector with the lower bounds (opt, default
%                               -inf for reversible reactions and 0 for
%                               irreversible)
%            ub                 vector with the upper bounds (opt, default
%                               inf)
%            c                  vector with the objective function
%                               coefficients (opt, default 0)
%            eccodes            cell array with the EC-numbers for each
%                               reactions. Delimit several EC-numbers with
%                               ";" (opt, default '')
%            subSystems         cell array with the subsystems for each
%                               reaction (opt, default '')
%            grRules            cell array with the gene-reaction
%                               relationship for each reaction. For example
%                               "(A and B) or (C)" means that the reaction
%                               could be catalyzed by a complex between
%                               A & B or by C on its own. All the genes
%                               have to be present in model.genes. Add
%                               genes with addGenesRaven before calling
%                               this function if needed (opt, default '')
%            rxnMiriams         cell array with Miriam structures (opt,
%                               default [])
%            rxnComps           cell array with compartments (as in
%                               model.comps) (opt, default {})
%            rxnNotes           cell array with reaction notes (opt,
%                               default '')
%            rxnReferences      cell array with reaction references (opt,
%                               default '')
%            rxnConfidenceScores   vector with reaction confidence scores
%                               (opt, default NaN)
%   eqnType          double describing how the equation string should be
%                    interpreted
%                    1 - The metabolites are matched to model.mets. New
%                        metabolites (if allowed) are added to
%                        "compartment"
%                    2 - The metabolites are matched to model.metNames and
%                        all metabolites are assigned to "compartment". Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are
%                        already used in the model then the numbering will
%                        start from the highest used integer+1
%                    3 - The metabolites are written as
%                        "metNames[comps]". Only compartments in
%                        model.comps are allowed. Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are
%                        already used in the model then the numbering will
%                        start from the highest used integer+1
%   compartment      a string with the compartment the metabolites should
%                    be placed in when using eqnType=2. Must match
%                    model.comps (opt when eqnType=1 or eqnType=3)
%   allowNewMets     true if the function is allowed to add new
%                    metabolites. It is highly recommended to first add
%                    any new metabolites with addMets rather than
%                    automatically through this function. addMets supports
%                    more annotation of metabolites, allows for the use of
%                    exchange metabolites, and using it reduces the risk
%                    of parsing errors (opt, default false)
%
%   newModel         an updated model structure
%
%   NOTE: This function does not make extensive checks about formatting of
%   gene-reaction rules.
%
%   NOTE: When adding metabolites to a compartment where it previously
%   doesn't exist, the function will copy any available information from
%   the metabolite in another compartment.
%
%   Usage: newModel=addRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
%
%   Simonas Marcisauskas, 2018-04-03
%

if nargin<4
    compartment=[];
end
if nargin<5
    allowNewMets=false;
end

newModel=model;

%If no reactions should be added
if isempty(rxnsToAdd)
    return;
end

%Check the input
if ~isnumeric(eqnType)
    EM='eqnType must be numeric';
    dispEM(EM);
else
    if ~ismember(eqnType,[1 2 3])
        EM='eqnType must be 1, 2, or 3';
        dispEM(EM);
    end
end

if eqnType==2 || (eqnType==1 && allowNewMets==true)
    if ~ischar(compartment)
        EM='compartment must be a string';
        dispEM(EM);
    end
    if ~ismember(compartment,model.comps)
        EM='compartment must match one of the compartments in model.comps';
        dispEM(EM);
    end
end

if ~isfield(rxnsToAdd,'rxns')
    EM='rxns is a required field in rxnsToAdd';
    dispEM(EM);
elseif iscell(rxnsToAdd.rxns)
    %To fit with some later printing
    rxnsToAdd.rxns=rxnsToAdd.rxns(:);
end
if ~isfield(rxnsToAdd,'equations')
    EM='equations is a required field in rxnsToAdd';
    dispEM(EM);
end

if any(ismember(rxnsToAdd.rxns,model.rxns))
    EM='One or more reaction id was already present in the model. Use changeRxns or remove the existing reactions before adding new ones';
    dispEM(EM);
end

if ~iscellstr(rxnsToAdd.rxns) && ~ischar(rxnsToAdd.rxns)
    %It could also be a string, but it's not encouraged
    EM='rxnsToAdd.rxns must be a cell array of strings';
    dispEM(EM);
else
    rxnsToAdd.rxns=cellstr(rxnsToAdd.rxns);
end
if ~iscellstr(rxnsToAdd.equations) && ~ischar(rxnsToAdd.equations)
    %It could also be a string, but it's not encouraged
    EM='rxnsToAdd.equations must be a cell array of strings';
    dispEM(EM);
else
    rxnsToAdd.equations=cellstr(rxnsToAdd.equations);
end

nRxns=numel(rxnsToAdd.rxns);
nOldRxns=numel(model.rxns);
filler=cell(nRxns,1);
filler(:)={''};
largeFiller=cell(nOldRxns,1);
largeFiller(:)={''};

%***Add everything to the model except for the equations.
if numel(rxnsToAdd.equations)~=nRxns
    EM='rxnsToAdd.equations must have the same number of elements as rxnsToAdd.rxns';
    dispEM(EM);
end

%Parse the equations. This is done at this early stage since I need the
%reversibility info
[S, mets, badRxns, reversible]=constructS(rxnsToAdd.equations);
EM='The following equations have one or more metabolites both as substrate and product. Only the net equations will be added:';
dispEM(EM,false,rxnsToAdd.rxns(badRxns));

newModel.rev=[newModel.rev;reversible];
newModel.rxns=[newModel.rxns;rxnsToAdd.rxns(:)];

if isfield(rxnsToAdd,'rxnNames')
    if numel(rxnsToAdd.rxnNames)~=nRxns
        EM='rxnsToAdd.rxnNames must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnNames')
        newModel.rxnNames=largeFiller;
    end
    newModel.rxnNames=[newModel.rxnNames;rxnsToAdd.rxnNames(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnNames')
        newModel.rxnNames=[newModel.rxnNames;filler];
    end
end

if isfield(rxnsToAdd,'lb')
    if numel(rxnsToAdd.lb)~=nRxns
        EM='rxnsToAdd.lb must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'lb')
        newModel.lb=zeros(nOldRxns,1);
        newModel.lb(newModel.rev~=0)=-inf;
    end
    newModel.lb=[newModel.lb;rxnsToAdd.lb(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'lb')
        I=zeros(nRxns,1);
        I(reversible~=0)=-inf;
        newModel.lb=[newModel.lb;I];
    end
end

if isfield(rxnsToAdd,'ub')
    if numel(rxnsToAdd.ub)~=nRxns
        EM='rxnsToAdd.ub must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'ub')
        newModel.ub=inf(nOldRxns,1);
    end
    newModel.ub=[newModel.ub;rxnsToAdd.ub(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'ub')
        newModel.ub=[newModel.ub;inf(nRxns,1)];
    end
end

if isfield(rxnsToAdd,'c')
    if numel(rxnsToAdd.c)~=nRxns
        EM='rxnsToAdd.c must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'c')
        newModel.c=zeros(nOldRxns,1);
    end
    newModel.c=[newModel.c;rxnsToAdd.c(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'c')
        newModel.c=[newModel.c;zeros(nRxns,1)];
    end
end

if isfield(rxnsToAdd,'eccodes')
    if numel(rxnsToAdd.eccodes)~=nRxns
        EM='rxnsToAdd.eccodes must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'eccodes')
        newModel.eccodes=largeFiller;
    end
    newModel.eccodes=[newModel.eccodes;rxnsToAdd.eccodes(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'eccodes')
        newModel.eccodes=[newModel.eccodes;filler];
    end
end

if isfield(rxnsToAdd,'subSystems')
    if numel(rxnsToAdd.subSystems)~=nRxns
        EM='rxnsToAdd.subSystems must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'subSystems')
        newModel.subSystems=largeFiller;
    end
    newModel.subSystems=[newModel.subSystems;rxnsToAdd.subSystems(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'subSystems')
        newModel.subSystems=[newModel.subSystems;filler];
    end
end
if isfield(rxnsToAdd,'rxnMiriams')
    if numel(rxnsToAdd.rxnMiriams)~=nRxns
        EM='rxnsToAdd.rxnMiriams must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnMiriams')
        newModel.rxnMiriams=cell(nOldRxns,1);
    end
    newModel.rxnMiriams=[newModel.rxnMiriams;rxnsToAdd.rxnMiriams(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnMiriams')
        newModel.rxnMiriams=[newModel.rxnMiriams;cell(nRxns,1)];
    end
end
if isfield(rxnsToAdd,'rxnComps')
    if numel(rxnsToAdd.rxnComps)~=nRxns
        EM='rxnsToAdd.rxnComps must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnComps')
        newModel.rxnComps=ones(nOldRxns,1);
        EM='Adding reactions with compartment information to a model without such information. All existing reactions will be assigned to the first compartment';
        dispEM(EM,false);
    end
    newModel.rxnComps=[newModel.rxnComps;rxnsToAdd.rxnComps(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnComps')
        newModel.rxnComps=[newModel.rxnComps;ones(nRxns,1)];
        %fprintf('NOTE: The added reactions will be assigned to the first
        %compartment\n');
    end
end
if isfield(rxnsToAdd,'grRules')
    if numel(rxnsToAdd.grRules)~=nRxns
        EM='rxnsToAdd.grRules must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'grRules')
        newModel.grRules=largeFiller;
    end
    newModel.grRules=[newModel.grRules;rxnsToAdd.grRules(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'grRules')
        newModel.grRules=[newModel.grRules;filler];
    end
end

if isfield(rxnsToAdd,'rxnFrom')
    if numel(rxnsToAdd.rxnFrom)~=nRxns
        EM='rxnsToAdd.rxnFrom must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnFrom')
        newModel.rxnFrom=largeFiller;
    end
    newModel.rxnFrom=[newModel.rxnFrom;rxnsToAdd.rxnFrom(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnFrom')
        newModel.rxnFrom=[newModel.rxnFrom;filler];
    end
end

if isfield(rxnsToAdd,'rxnNotes')
    if numel(rxnsToAdd.rxnNotes)~=nRxns
        EM='rxnsToAdd.rxnNotes must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnNotes')
        newModel.rxnNotes=largeFiller;
    end
    newModel.rxnNotes=[newModel.rxnNotes;rxnsToAdd.rxnNotes(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnNotes')
        newModel.rxnNotes=[newModel.rxnNotes;filler];
    end
end

if isfield(rxnsToAdd,'rxnReferences')
    if numel(rxnsToAdd.rxnReferences)~=nRxns
        EM='rxnsToAdd.rxnReferences must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnReferences')
        newModel.rxnReferences=largeFiller;
    end
    newModel.rxnReferences=[newModel.rxnReferences;rxnsToAdd.rxnReferences(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnReferences')
        newModel.rxnReferences=[newModel.rxnReferences;filler];
    end
end

if isfield(rxnsToAdd,'pwys')
    if numel(rxnsToAdd.pwys)~=nRxns
        EM='rxnsToAdd.pwys must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'pwys')
        newModel.pwys=largeFiller;
    end
    newModel.pwys=[newModel.pwys;rxnsToAdd.pwys(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'pwys')
        newModel.pwys=[newModel.pwys;filler];
    end
end

if isfield(rxnsToAdd,'rxnConfidenceScores')
    if numel(rxnsToAdd.rxnConfidenceScores)~=nRxns
        EM='rxnsToAdd.rxnConfidenceScores must have the same number of elements as rxnsToAdd.rxns';
        dispEM(EM);
    end
    %Fill with standard if it doesn't exist
    if ~isfield(newModel,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=NaN(nOldRxns,1);
        EM='Adding reactions with confidence scores without such information. All existing reactions will have confidence scores as NaNs';
        dispEM(EM,false);
    end
    newModel.rxnConfidenceScores=[newModel.rxnConfidenceScores;rxnsToAdd.rxnConfidenceScores(:)];
else
    %Fill with standard if it doesn't exist
    if isfield(newModel,'rxnConfidenceScores')
        newModel.rxnConfidenceScores=[newModel.rxnConfidenceScores;NaN(nRxns,1)];
    end
end


%***Start parsing the equations and adding the info to the S matrix The
%mets are matched to model.mets
if eqnType==1
    [I, J]=ismember(mets,model.mets);
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.mets=mets(~I);
            metsToAdd.metNames=metsToAdd.mets;
            metsToAdd.compartments=compartment;
            newModel=addMets(newModel,metsToAdd);
        else
            EM='One or more equations contain metabolites that are not in model.mets. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function';
            dispEM(EM);
        end
    end
    %Calculate the indexes of the metabolites and add the info
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%Do some stuff that is the same for eqnType=2 and eqnType=3
if eqnType==2 || eqnType==3
    %For later..
    t2=strcat(model.metNames,'***',model.comps(model.metComps));
end

%The mets are matched to model.metNames and assigned to "compartment"
if eqnType==2
    %%Check that the metabolite names aren't present in the same
    %%compartment.
    %Not the neatest way maybe..
    t1=strcat(mets,'***',compartment);
    [I, J]=ismember(t1,t2);
    
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.metNames=mets(~I);
            metsToAdd.compartments=compartment;
            newModel=addMets(newModel,metsToAdd);
        else
            EM='One or more equations contain metabolites that are not in model.mets. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function';
            dispEM(EM);
        end
    end
    
    %Calculate the indexes of the metabolites
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%The equations are on the form metNames[compName]
if eqnType==3
    %Parse the metabolite names
    metNames=cell(numel(mets),1);
    compartments=metNames;
    for i=1:numel(mets)
        starts=max(strfind(mets{i},'['));
        ends=max(strfind(mets{i},']'));
        
        %Check that the formatting is correct
        if isempty(starts) || isempty(ends) || ends<numel(mets{i})
            EM=['The metabolite ' mets{i} ' is not correctly formatted for eqnType=3'];
            dispEM(EM);
        end
        
        %Check that the compartment is correct
        compartments{i}=mets{i}(starts+1:ends-1);
        I=ismember(compartments(i),newModel.comps);
        if ~I
            EM=['The metabolite ' mets{i} ' has a compartment that is not in model.comps'];
            dispEM(EM);
        end
        metNames{i}=mets{i}(1:starts-1);
    end
    
    %Check if the metabolite exists already
    t1=strcat(metNames,'***',compartments);
    [I, J]=ismember(t1,t2);
    
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.metNames=metNames(~I);
            metsToAdd.compartments=compartments(~I);
            newModel=addMets(newModel,metsToAdd);
        else
            EM='One or more equations contain metabolites that are not in model.metNames. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function';
            dispEM(EM);
        end
    end
    
    %Calculate the indexes of the metabolites
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%Add the info to the stoichiometric matrix
newModel.S=[newModel.S sparse(size(newModel.S,1),nRxns)];
for i=1:nRxns
    newModel.S(metIndexes,nOldRxns+i)=S(:,i);
    %Parse the grRules and check whether all genes in grRules appear in
    %genes
    if isfield(newModel,'grRules')
        rule=newModel.grRules{nOldRxns+i};
        rule=strrep(rule,'(','');
        rule=strrep(rule,')','');
        rule=strrep(rule,' or ',' ');
        rule=strrep(rule,' and ',' ');
        genes=regexp(rule,' ','split');
        [I, J]=ismember(genes,newModel.genes);
        if ~all(I) && any(rule)
            EM=['Not all genes for reaction ' rxnsToAdd.rxns{i} ' were found in model.genes. If needed, add genes with addGenesRaven before calling this function'];
            dispEM(EM);
        end
    end
end

%Fix grRules and reconstruct rxnGeneMat
[grRules,rxnGeneMat] = standardizeGrRules(newModel,true);
newModel.grRules = grRules;
newModel.rxnGeneMat = rxnGeneMat;
end
