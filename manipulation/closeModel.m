function closedModel=closeModel(model)
% closeModel  Add boundary metabolites and their exchange reactions.
%
% Adds boundary metabolites and their participation in exchange reactions.
%
% Parameters
% ----------
% model : struct
%     a model structure.
%
% Returns
% -------
% closedModel : struct
%     an updated model structure with boundary metabolites added.
%
% Examples
% --------
%     closedModel = closeModel(model);

closedModel=model;

closedModel.comps{numel(closedModel.comps)+1,1}='b';
closedModel.compNames{numel(closedModel.compNames)+1,1}='boundary';
if isfield(closedModel,'compMiriams')
    closedModel.compMiriams{numel(closedModel.compMiriams)+1,1}=[];
end
closedModel.unconstrained=zeros(numel(closedModel.mets),1);

for i=1:numel(closedModel.rxns)
    col=closedModel.S(:,i);
    %A boundary/exchange/sink/demand reaction has metabolites on only one
    %side, same definition as getExchangeRxns: no substrates or no
    %products. The previous coefficient-magnitude rule (all nonzero
    %coefficients summing to 1 in absolute value) missed a scaled
    %single-metabolite sink such as "2 A =>" (sum 2) outright, and only
    %matched a multi-metabolite one-sided reaction such as
    %"0.5 A + 0.5 B =>" by coincidence, when its magnitudes summed to
    %exactly 1 --- and even then mishandled it below, since find() on such
    %a reaction returns more than one index and only the first metabolite
    %was ever used to build the boundary metabolite, silently dropping the
    %rest.
    if any(col>0) && any(col<0)
        continue
    end
    metIdx=find(col);
    if isempty(metIdx)
        continue
    end
    %A single metabolite names and inherits properties from that
    %metabolite, as before. Several metabolites on the one populated side
    %have no single metabolite to name the boundary after or copy
    %properties from, so it is named after the reaction instead and left
    %without those properties.
    oneMet=isscalar(metIdx);
    if oneMet
        newMetId=strcat(closedModel.mets{metIdx},'_b');
    else
        newMetId=strcat(closedModel.rxns{i},'_b');
    end
    closedModel.mets{numel(closedModel.mets)+1}=newMetId;
    if isfield(closedModel,'metNames')
        if oneMet
            closedModel.metNames{numel(closedModel.metNames)+1}=closedModel.metNames{metIdx};
        else
            closedModel.metNames{numel(closedModel.metNames)+1}=newMetId;
        end
    end
    closedModel.metComps(numel(closedModel.metComps)+1)=numel(closedModel.comps);
    if isfield(closedModel,'inchis')
        if oneMet
            closedModel.inchis{numel(closedModel.inchis)+1}=closedModel.inchis{metIdx};
        else
            closedModel.inchis{numel(closedModel.inchis)+1}='';
        end
    end
    if isfield(closedModel,'metSmiles')
        if oneMet
            closedModel.metSmiles{numel(closedModel.metSmiles)+1}=closedModel.metSmiles{metIdx};
        else
            closedModel.metSmiles{numel(closedModel.metSmiles)+1}='';
        end
    end
    if isfield(closedModel,'metFormulas')
        if oneMet
            closedModel.metFormulas{numel(closedModel.metFormulas)+1}=closedModel.metFormulas{metIdx};
        else
            closedModel.metFormulas{numel(closedModel.metFormulas)+1}='';
        end
    end
    if isfield(closedModel,'metMiriams')
        if oneMet
            closedModel.metMiriams{numel(closedModel.metMiriams)+1}=closedModel.metMiriams{metIdx};
        else
            closedModel.metMiriams{numel(closedModel.metMiriams)+1}=[];
        end
    end
    if isfield(closedModel,'metFrom')
        closedModel.metFrom{numel(closedModel.metFrom)+1}='';
    end
    if isfield(closedModel,'metCharges')
        if oneMet
            closedModel.metCharges(numel(closedModel.metCharges)+1)=closedModel.metCharges(metIdx);
        else
            closedModel.metCharges(numel(closedModel.metCharges)+1)=0;
        end
    end
    if isfield(closedModel,'metDeltaG')
        if oneMet
            closedModel.metDeltaG(numel(closedModel.metDeltaG)+1)=closedModel.metDeltaG(metIdx);
        else
            closedModel.metDeltaG(numel(closedModel.metDeltaG)+1)=NaN;
        end
    end
    closedModel.unconstrained(numel(closedModel.unconstrained)+1)=1;
    closedModel.b(numel(closedModel.b)+1)=0;
    closedModel.S=[closedModel.S;sparse(1,size(closedModel.S,2))];
    %Any nonzero value blocks the reaction once its slack is pinned to 0 by
    %model.b above, regardless of magnitude; the sign is chosen only so the
    %boundary row reads the same way a human would expect (produced by vs.
    %consumed by the reaction).
    if sum(col)>0
        closedModel.S(numel(closedModel.mets),i)=-1;
    else
        closedModel.S(numel(closedModel.mets),i)=1;
    end
end

end
