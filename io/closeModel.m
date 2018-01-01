function closedModel=closeModel(model)
% closeModel
%   Adds boundary metabolites and their participation in exchange
%   reactions.
%
%   model                 a model structure
%
%   closedModel           an updated closedModel structure
%
%   Usage: closedModel=closeModel(model)
%
%   Simonas Marcisauskas, 2017-09-18
%


closedModel=model;

closedModel.comps{numel(closedModel.comps)+1}='b';
closedModel.compNames{numel(closedModel.compNames)+1}='boundary';
if isfield(closedModel,'compMiriams')
    closedModel.compMiriams{numel(closedModel.compMiriams)+1}=[];
end;
closedModel.unconstrained=zeros(numel(closedModel.mets),1);

for i=1:numel(closedModel.rxns)
    if sum(abs(closedModel.S(:,i)))==1
        closedModel.mets{numel(closedModel.mets)+1}=strcat(closedModel.mets{find(closedModel.S(:,i))},'_b');
        if isfield(closedModel,'metNames')
            closedModel.metNames{numel(closedModel.metNames)+1}=closedModel.metNames{find(closedModel.S(:,i))};
        end;
        closedModel.metComps(numel(closedModel.metComps)+1)=numel(closedModel.comps);
        if isfield(closedModel,'inchis')
            closedModel.inchis{numel(closedModel.inchis)+1}=closedModel.inchis{find(closedModel.S(:,i))};
        end;
        if isfield(closedModel,'metFormulas')
            closedModel.metFormulas{numel(closedModel.metFormulas)+1}=closedModel.metFormulas{find(closedModel.S(:,i))};
        end;
        if isfield(closedModel,'metMiriams')
            closedModel.metMiriams{numel(closedModel.metMiriams)+1}=closedModel.metMiriams{find(closedModel.S(:,i))};
        end;
        if isfield(closedModel,'metCharges')
            closedModel.metCharges(numel(closedModel.metCharges)+1)=closedModel.metCharges(find(closedModel.S(:,i)));
        end;
        closedModel.unconstrained(numel(closedModel.unconstrained)+1)=1;
        closedModel.b(numel(closedModel.b)+1)=0;
        closedModel.S=[closedModel.S;sparse(1,size(closedModel.S,2))];
        if sum(closedModel.S(:,i))==1
            closedModel.S(numel(closedModel.mets),i)=-1;
        else
            closedModel.S(numel(closedModel.mets),i)=1;
        end;
    end;
end;

end
