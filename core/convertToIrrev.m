function irrevModel=convertToIrrev(model,rxns)
% convertToIrrev
%   Converts a model to irreversible form
%
%   model         a model structure
%   rxns          cell array with the reactions so split (if reversible)
%                 (opt, default model.rxns)
%
%   irrevModel    a model structure where reversible reactions have
%                 been split into one forward and one reverse reaction
%
%   The reverse reactions are saved as 'rxnID_REV'.
%
%   Usage: irrevModel=convertToIrrev(model,rxns)
%
%   Rasmus Agren, 2013-08-01
%

if nargin<2
    rxns=model.rxns;
end

irrevModel=model;

I=getIndexes(model,rxns,'rxns',true);

revIndexesBool=model.rev~=0 & I;
revIndexes=find(revIndexesBool);
if any(revIndexesBool)
    irrevModel.S=[model.S,model.S(:,revIndexes)*-1];
    irrevModel.rev(revIndexes)=0;
    irrevModel.rev=[irrevModel.rev;zeros(numel(revIndexes),1)];

    %Get the limits for all normal and reversible rxns
    ubNormal=irrevModel.ub;
    ubNormal(revIndexes(ubNormal(revIndexes)<0))=0;
    lbNormal=irrevModel.lb;
    lbNormal(revIndexes(lbNormal(revIndexes)<0))=0;
    ubRev=irrevModel.lb(revIndexes)*-1;
    ubRev(ubRev<0)=0;
    lbRev=irrevModel.ub(revIndexes)*-1;
    lbRev(lbRev<0)=0;
    irrevModel.ub=[ubNormal;ubRev];
    irrevModel.lb=[lbNormal;lbRev];

    %The objective coefficents should be zero for the backwards reversible
    %reactions unless they were negative in the original. In that case they
    %should be positive for the backwards reversible and deleted from the
    %original
    irrevC=zeros(numel(revIndexes),1);

    if any(irrevModel.c(revIndexes)<0)
        originalC=irrevModel.c(revIndexes);
        irrevC(irrevModel.c(revIndexes)<0)=originalC(originalC<0)*-1;
        irrevModel.c(irrevModel.c<0 & revIndexesBool)=0;
    end
    irrevModel.c=[irrevModel.c;irrevC];

    irrevModel.rxns=[irrevModel.rxns;strcat(irrevModel.rxns(revIndexes),'_REV')];
    irrevModel.rxnNames=[irrevModel.rxnNames;strcat(irrevModel.rxnNames(revIndexes),' (reversible)')];

    if isfield(irrevModel,'grRules')
        irrevModel.grRules=[irrevModel.grRules;irrevModel.grRules(revIndexes,:)];
    end
    if isfield(irrevModel,'rxnMiriams')
        irrevModel.rxnMiriams=[irrevModel.rxnMiriams;irrevModel.rxnMiriams(revIndexes,:)];
    end
    if isfield(irrevModel,'rxnGeneMat')
        irrevModel.rxnGeneMat=[irrevModel.rxnGeneMat;irrevModel.rxnGeneMat(revIndexes,:)];
    end
    if isfield(irrevModel,'subSystems')
    	irrevModel.subSystems=[irrevModel.subSystems;irrevModel.subSystems(revIndexes)];
    end
    if isfield(irrevModel,'eccodes')
    	irrevModel.eccodes=[irrevModel.eccodes;irrevModel.eccodes(revIndexes)];
    end
    if isfield(irrevModel,'rxnComps')
    	irrevModel.rxnComps=[irrevModel.rxnComps;irrevModel.rxnComps(revIndexes)];
    end
    if isfield(irrevModel,'rxnFrom')
    	irrevModel.rxnFrom=[irrevModel.rxnFrom;irrevModel.rxnFrom(revIndexes)];
    end
    if isfield(irrevModel,'rxnScores')
    	irrevModel.rxnScores=[irrevModel.rxnScores;irrevModel.rxnScores(revIndexes)];
    end
    if isfield(irrevModel,'rxnNotes')
    	irrevModel.rxnNotes=[irrevModel.rxnNotes;irrevModel.rxnNotes(revIndexes)];
    end
    if isfield(irrevModel,'rxnConfidenceScores')
    	irrevModel.rxnConfidenceScores=[irrevModel.rxnConfidenceScores;irrevModel.rxnConfidenceScores(revIndexes)];
    end
end
end
