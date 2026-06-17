function targets=FSEOF(model,biomassRxn,targetRxn,varargin)
% FSEOF  Flux Scanning based on Enforced Objective Flux.
%
% Implements the Flux Scanning based on Enforced Objective Flux algorithm.
% This function writes a tab-delimited file or prints to the command
% window. If an output has been specified (targets), it will also generate
% a structure indicating for each model reaction whether it is identified
% by FSEOF as a target and the slope of the reaction when switching from
% biomass formation to product formation.
%
% Parameters
% ----------
% model : struct
%     a model structure.
% biomassRxn : char
%     reaction ID of the biomass formation or growth reaction.
% targetRxn : char
%     reaction ID of the target reaction.
%
% Name-Value Arguments
% --------------------
% iterations : double
%     number of iterations (default 10).
% coefficient : double
%     ratio of optimal target reaction flux, must be less than 1
%     (default 0.9).
% outputFile : char
%     output filename (default prints to command window).
% corrThreshold : double
%     minimum absolute Pearson correlation coefficient between a
%     reaction's absolute flux and the enforced product flux for the
%     reaction to be selected as a target. Values in [0, 1]; higher
%     values require a more linear response and suppress LP alternative-
%     optima noise (default 0.9).
%
% Returns
% -------
% targets : struct
%     structure with information for identified targets, with fields:
%
%     - logical : logical array, true for amplification targets (backward
%       compatible alias for amplify).
%     - amplify : logical array, reactions whose absolute flux increases
%       with enforced product flux.
%     - knockdown : logical array, reactions driven toward but not to zero.
%     - knockout : logical array, reactions driven to zero at max product.
%     - slope : numeric array with per-reaction regression slopes of
%       absolute flux versus enforced product flux (all reactions, not
%       just targets).
%
% Examples
% --------
%     targets = FSEOF(model, biomassRxn, targetRxn, iterations, ...
%                     coefficient, outputFile);

biomassRxn=char(biomassRxn);
targetRxn=char(targetRxn);

p=parseRAVENargs(varargin, {'iterations',10; 'coefficient',0.9; 'outputFile',[]; 'corrThreshold',0.9});
iterations=p.iterations;
coefficient=p.coefficient;
outputFile=p.outputFile;
corrThreshold=p.corrThreshold;

%Find out the maximum theoretical yield of target reaction
model=setParam(model,'obj',targetRxn,1);
sol=solveLP(model,1);
targetMax=sol.f*coefficient;

model=setParam(model,'obj',biomassRxn,1);

nRxns=length(model.rxns);
results=zeros(nRxns,iterations);
enforcedFlux=(1:iterations)*targetMax/iterations;

%Run pFBA at each enforced target flux level
for i=1:iterations
    model=setParam(model,'lb',targetRxn,enforcedFlux(i));
    sol=solveLP(model,1);
    results(:,i)=sol.x;
end

%Compute per-reaction regression slope (FS4) and Pearson correlation (FS1)
%on absolute flux values to treat forward/reverse symmetrically.
absResults=abs(results);
slope=zeros(nRxns,1);
rCoeff=zeros(nRxns,1);
hasVariation = iterations > 1 && std(enforcedFlux) > 0;
if hasVariation
    for num=1:nRxns
        fluxTrack=absResults(num,:);
        if any(fluxTrack > 1e-10) && std(fluxTrack) > 0
            p_fit=polyfit(enforcedFlux,fluxTrack,1);
            slope(num)=p_fit(1);
            r=corrcoef(enforcedFlux,fluxTrack);
            rCoeff(num)=r(1,2);
        end
    end
end

%Classify targets: select reactions whose absolute flux is sufficiently
%correlated with the enforced product flux (FS1), then split by direction
%and whether any flux remains at maximum enforcement (FS2).
isCorrelated = abs(rCoeff) >= corrThreshold;
isAmplify    = isCorrelated & slope > 0;
isKnockdown  = isCorrelated & slope < 0 & absResults(:,end) > 1e-8;
isKnockout   = isCorrelated & slope < 0 & absResults(:,end) <= 1e-8;

rxnDirection=sign(results(:,1));

%Generating output
formatSpec='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
if ~isempty(outputFile)
    outputFile=char(outputFile);
    fid=fopen(outputFile,'w');
    fprintf(fid,formatSpec,'Type','Slope','rowID','Enzyme ID','Enzyme Name','Subsystems','Direction','Gr Rule');
else
    fid=-1;
    fprintf(formatSpec,'Type','Slope','rowID','Enzyme ID','Enzyme Name','Subsystems','Direction','Gr Rule');
end

targetTypes  = {'Amplify','Knockdown','Knockout'};
targetArrays = {isAmplify, isKnockdown, isKnockout};

for t=1:3
    flag=targetArrays{t};
    for num=1:nRxns
        if flag(num)
            A_type=targetTypes{t};
            A0=num2str(slope(num));
            A1=num2str(num);
            A2=char(model.rxns(num));
            A3=char(model.rxnNames(num));
            if isfield(model,'subSystems') && ~isempty(model.subSystems{num})
                if ~any(cellfun(@(x) iscell(x), model.subSystems))
                    subSys=cellfun(@(x) {x}, model.subSystems, 'uni', 0);
                else
                    subSys=model.subSystems;
                end
                A4=char(strjoin(subSys{num},';'));
            else
                A4='';
            end
            A5=num2str(model.rev(num)*rxnDirection(num));
            A6=char(model.grRules(num));
            if fid == -1
                fprintf(formatSpec,A_type,A0,A1,A2,A3,A4,A5,A6);
            else
                fprintf(fid,formatSpec,A_type,A0,A1,A2,A3,A4,A5,A6);
            end
        end
    end
end

if fid ~= -1
    fclose(fid);
end

if nargout == 1
    targets.logical   = isAmplify;
    targets.amplify   = isAmplify;
    targets.knockdown = isKnockdown;
    targets.knockout  = isKnockout;
    targets.slope     = slope;
end
end
