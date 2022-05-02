function [deletedRxns,metProduction,res,turnedOnRxns,fluxes]=ftINITInternalAlg(model,rxnScores,metData,essentialRxns,prodWeight,allowExcretion,remPosRev,params)
% ftINITInternalAlg
%	Generates a model using the INIT algorithm, based on proteomics and/or
%   transcriptomics and/or metabolomics and/or metabolic tasks
%
%   model           a reference model structure
%   rxnScores       a vector of scores for the reactions in the model.
%                   Positive scores are reactions to keep and negative
%                   scores are reactions to exclude (opt, default all 0.0)
%   metData         boolean matrix with mets as rows and rxns as columns
%                   saying which reaction produces each detected met (opt, default [])
%   essentialRxns   cell array of reactions that are essential and that
%                   have to be in the resulting model. This is normally
%                   used when fitting a model to task (see fitTasks) (opt,
%                   default [])
%   prodWeight      a score that determines the value of having
%                   net-production of metabolites. This is a way of having
%                   a more functional network as it provides a reason for
%                   including bad reactions for connectivity reasons. This
%                   score is for each metabolite, and the sum of these weights
%                   and the scores for the reactions is what is optimized
%                   (opt, default 0.5)
%   allowExcretion  true if excretion of all metabolites should be allowed.
%                   This results in fewer reactions being considered
%                   dead-ends, but all reactions in the resulting model may
%                   not be able to carry flux. If this is "false" then the
%                   equality constraints are taken from model.b. If the
%                   input model lacks exchange reactions then this should
%                   probably be "true", or a large proportion of the model
%                   would be excluded for connectivity reasons
%                   (opt, default false)
%   remPosRev       If true, the positive reversible reactions are removed from the problem.
%                   This is used in step 1 of ftINIT (opt, default false)
%   params          parameter structure as used by getMILPParams (opt,
%                   default [])
%
%   deletedRxns     reactions which were deleted by the algorithm
%   metProduction   array that indicates which of the
%                   metabolites in presentMets that could be
%                   produced
%                   -2: metabolite name not found in model
%                   -1: metabolite found, but it could not be produced
%                   1: metabolite could be produced
%   res             The result from the MILP
%   turnedOnRxns    The reactions determined to be present
%   fluxes          The fluxes from the MILP
%
%   This function is the actual implementation of the algorithm. See
%   getINITModel9 for a higher-level function for model reconstruction. 
%
%   Usage: [deletedRxns,metProduction,res,turnedOnRxns,fluxes]=runINIT9(model,...
%           rxnScores,presentMets,essentialRxns,prodWeight,allowExcretion,...
%           remPosRev,params)

if nargin<2
    rxnScores=zeros(numel(model.rxns),1);
end
if isempty(rxnScores)
    rxnScores=zeros(numel(model.rxns),1);
end
if nargin<3
    metData={};
end
if nargin<4
    essentialRxns={};
end
if isempty(essentialRxns)
    essentialRxns={};
end
essentialRxns=essentialRxns(:);
if nargin<5
    prodWeight=0.5;
end
if isempty(prodWeight)
    prodWeight=0.5;
end
if nargin<6
    allowExcretion = false;
end

if nargin<7
    remPosRev = false;
end

if nargin<8
    params=[];
end

%The model should be in the reversible format and all relevant exchange
%reactions should be open
if isfield(model,'unconstrained')
    EM='Exchange metabolites are still present in the model. Use simplifyModel if this is not intended';
    dispEM(EM,false);
end


%Get the indexes of the essential reactions and remove them from the
%scoring vector'
essential = ismember(model.rxns,essentialRxns);

essentialIndex=find(essential);


%Some nice to have numbers
nMets=numel(model.mets);
nRxns=numel(model.rxns);
nRxnsWithOnOff = nRxns;

%Reactions with score 0 will just be left in the model, and will not be part of the problem (but can carry flux).
%It is possible to set score = 0 for e.g. spontaneous reactions, exchange rxns, etc., which may not be that interesting
%to remove.

%if makeIrrev is on, we can just as well skip all positive reversible rxns - they will be turned on without carrying any flux
%since they can form a loop within themself (fwd-rev)
if remPosRev
    rxnScores(rxnScores > 0 & model.rev~=0) = 0;
end

%Handle metabolomics
%A problem with the metabolomics is that some of the producer reactions
%for a metabolite could be excluded from the problem. We solve this by
%adding them with score 0. They can still be seen as positive, since they
%either don't matter (there is another producer) or they can contribute to an 
%increased score.
revRxns=model.rev~=0; 
essRevRxns = find(revRxns & essential);
essIrrevRxns = find(~revRxns & essential);

if ~isempty(metData)
    %remove any metData rows that is connected to an essential rxn - 
    %these will be on regardless and will cause problems below.
    containsEssential = any(metData(:,essential),2);
    metData = metData(~containsEssential,:);
end

if ~isempty(metData)
    metRxns = any(metData,1).';
    posRxns = rxnScores > 0 | ((rxnScores == 0) & metRxns); 
else
    posRxns = rxnScores > 0;
end
negRxns = rxnScores < 0;

posRevRxns = find(posRxns & revRxns & ~essential);
negRevRxns = find(negRxns & revRxns & ~essential);
posIrrevRxns = find(posRxns & ~revRxns & ~essential);
negIrrevRxns = find(negRxns & ~revRxns & ~essential);
nPosRev = numel(posRevRxns);
nNegRev = numel(negRevRxns);
nPosIrrev = numel(posIrrevRxns);
nNegIrrev = numel(negIrrevRxns);
nEssRev = numel(essRevRxns);
nEssIrrev = numel(essIrrevRxns);
nMetabolMets = size(metData,1);

milpModel = model;

%These six categories need to be handled separately:
%
%PosIrrev (positive reaction score, irreversible): 
%flux >= 0.1*Yi, flux - 0.1*Yi - VPI == 0, 0 <= Yi(onoff) <= 1, 0 <= VPI <= Inf
%The nice thing with rxns with positive score is that they do not require a boolean. The optimizer will
%strive for maximizing the Yi here, so it can be a continuous variable, where we just force a flux if
%the variable is on.
%PosRev: 
% 1: Split up the flux into positive and negative: flux - vnrp + vnrn == 0, 0 <= vprp,vprn <= Inf.
% 2. Force one of the fluxes on: vprp + vprn >= 0.1*Yi, vprp + vprn - 0.1*Yi - vprv1 == 0. 0 <= Yi(onoff) <= 1, vprv1 >= 0
% 3. We need a bool to make sure that one of vprp and vprn are always zero:
%    vprb E {0,1}:  vprp <= 100*vprb (if bool is 0, vprp is zero): vprp - 100*vprb + vprv2 == 0, vprv2 >= 0
%    vprn <= (1-vprb)*100 (if bool is one, vprn is zero): vprn + 100*vprb + vprv3 == 0, -100 <= vprv3 <= inf
%Neg Irrev:
%For negative scores, we need to use an integer, I see no way around that.
% flux < 100*Yi, Yi E {0,1}. Flux - 100*Yi + vni == 0
%Neg Rev:
%To cases:
%A: Witout irrev model
% 1: Split up the flux into positive and negative: flux - vnrp + vnrn == 0, 0 <= vprp,vprn <= Inf.
% 2: Force the Yi (on/off) var on if the flux is on: vnrp + vnrn <= 0.1 * Yi, Yi E {0,1}: vnrp + vnrn - 0.1 * Yi + vnrv1 == 0, vnrv1 >= 0.
%B: Irrev model
% 1. We now have two reactions, but still just boolean variable. We know that
%    the flux can only be positive, so we just say that flux1 + flux2 <= 0.1 * Yi, Yi E {0,1}: flux1 + flux2 - 0.1 * Yi + vnrv1 == 0, vnrv1 >= 0.
%    We don't care if there is any loop - there is just no benefit for the objective to generate one, and it doesn't matter.
%Ess Irrev (essential rxn, i.e. forced to carry flux, irreversible): 
%flux >= 0.1 - solved with lb = 0.1.
%Ess Rev: 
% 1: Split up the flux into positive and negative: flux - verp + vern == 0, 0 <= verp,vern <= Inf.
% 2. Force one of the fluxes on: verp + vern >= 0.1, verp + vern - verv1 == 0. verv1 >= 0.1
% 3. We need a bool to make sure that one of verp and vern are always zero:
%    verb E {0,1}:  verp <= 100*verb (if bool is 0, verp is zero): verp - 100*verb + verv2 == 0, verv2 >= 0
%    vern <= (1-verb)*100 (if bool is one, vern is zero): vern + 100*verb + verv3 == 0, -100 <= verv3 <= inf

% Now metabolomics
% The basic idea is that the bonus is gotten if any of the producer reactions are on. We can
% therefore use the "on" variables directly for the reactions that are included. We therefore add one variable per met, mon, which can be 
% continuous and between 0 and 1. We then say that mon <= v1 + v2 + ... + vn, i.e. that 
% mon + mv1 - v1 - v2 ... - vn == 0, 0 <= mon <= 1, 0 <= mv1 <= inf
% mon gets -prodWeight in the c vector (i.e. objective function)
% A tricky part is that some of the reactions producing a metabolite are left outside the problem
% This is solved by moving them into the problem, but with the score 0. They can still be treated as positive
% reactions. In the end, we are not interested if they are on though, but they may be
% allow for production of a metabolite, giving no benefit of turning on another producer reaction.

%The total matrix then looks like this. Note that we have ordered the categories, two reactions each,
%in the S matrix to make it easier to follow the figure. In practice, they come in random order, which
%is why the SEye variable is used further down.
%all categories below have two columns - the two column starts where label starts
%since the label cannot fit in two chars, they are on multiple lines
%
%The met variables and constraints are not in this figure for practical reasons.
%
%                           vprb  vprv3 vnrn  vern  verv2 mv1
% pi  ni  ei  Ypi Yni vpi vprn  vprv2 vnrp  verp  verb  mon
%   pr  nr  er  Ypr Ynr vprp  vprv1 vni   vnrv1 verv1 verv3
% SSSSSSSSSSSS0000000000000000000000000000000000000000000000
% SSSSSSSSSSSS0000000000000000000000000000000000000000000000 S block
% SSSSSSSSSSSS0000000000000000000000000000000000000000000000
% 1           N       -                                      Pos irrev block
%  1           N       -
%   1                   - 1                                  Pos rev block 1
%    1                   - 1
%               N       1 1   -                              Pos rev block 2
%                N       1 1   -
%                       1   M   1                            Pos rev block 3
%                        1   M   1
%                         1 C     1                          Pos rev block 4
%                          1 C     1
%     1           M                 1                        Neg irrev block
%      1           M                 1
%       1                             - 1                    Neg rev block 1
%        1                             - 1
%                   M                 1 1 1                  Neg rev block 2
%                    M                 1 1 1
%           1                               - 1              Ess rev block 1
%            1                               - 1
%                                           1 1 -            Ess rev block 2
%                                            1 1 -
%                                           1     M 1        Ess rev block 3
%                                            1     M 1
%                                             1   C   1      Ess rev block 4
%                                              1   C   1
%             --------                                  1 1  Met block - Here, we assume that all variables support each met.
%             --------                                   1 1             In practice, fewer of the "on" variables with -1 should be included
%
%
%M = -100
%N = -0.1
%D = 0.1
%C = 100
%- = -1
%
%It looks a bit different for the irrev model, but still quite similar

%build the A matrix
%S row
forceOnLim=0.0001;

%if makeIrrev
%    varsPerNegRev = 1;
%else
    varsPerNegRev = 3;
%end


%Figure out the number of variables needed for metabolomics
nMetVars = 2*size(metData,1); %mon and mv1

sRow = [milpModel.S sparse(nMets, nPosIrrev*2+nPosRev*7+nNegIrrev*2+nNegRev*(1+varsPerNegRev)+nEssRev*6 + nMetVars)];
sEye = speye(nRxns);
nYBlock = nPosIrrev+nPosRev+nNegIrrev+nNegRev;
piRows = [sEye(posIrrevRxns,:) speye(nPosIrrev)*-forceOnLim sparse(nPosIrrev,nPosRev+nNegIrrev+nNegRev) speye(nPosIrrev)*-1 sparse(nPosIrrev,nPosRev*6+nNegIrrev+nNegRev*varsPerNegRev+nEssRev*6+nMetVars)];

prRows1 = [sEye(posRevRxns,:) sparse(nPosRev,nYBlock + nPosIrrev) speye(nPosRev)*-1 speye(nPosRev) sparse(nPosRev,nPosRev*4+nNegIrrev+nNegRev*varsPerNegRev+nEssRev*6+nMetVars)];
prRows2 = [sparse(nPosRev,nRxns + nPosIrrev) speye(nPosRev)*-forceOnLim sparse(nPosRev,nNegIrrev + nNegRev + nPosIrrev) speye(nPosRev) speye(nPosRev) sparse(nPosRev,nPosRev) speye(nPosRev)*-1 sparse(nPosRev,nPosRev*2+nNegIrrev+nNegRev*varsPerNegRev+nEssRev*6+nMetVars)];
prRows3 = [sparse(nPosRev,nRxns + nYBlock + nPosIrrev) speye(nPosRev) sparse(nPosRev,nPosRev) speye(nPosRev)*-100 sparse(nPosRev,nPosRev) speye(nPosRev) sparse(nPosRev,nPosRev+nNegIrrev+nNegRev*varsPerNegRev+nEssRev*6+nMetVars)];
prRows4 = [sparse(nPosRev,nRxns + nYBlock + nPosIrrev + nPosRev) speye(nPosRev) speye(nPosRev)*100 sparse(nPosRev,nPosRev*2) speye(nPosRev) sparse(nPosRev, nNegIrrev+nNegRev*varsPerNegRev+nEssRev*6+nMetVars)];

niRows = [sEye(negIrrevRxns,:) sparse(nNegIrrev,nPosIrrev + nPosRev) speye(nNegIrrev)*-100 sparse(nNegIrrev,nNegRev + nPosIrrev + nPosRev*6) speye(nNegIrrev) sparse(nNegIrrev,nNegRev*varsPerNegRev + nEssRev*6+nMetVars)];

%if makeIrrev
    %flux1 + flux2 - 100 * Yi + vnrv1 == 0, vnrv1 >= 0.
%    nrRows = [sEye(negRevRxns,:) + sEye(reversedNegRevInd,:) sparse(nNegRev,nPosIrrev + nPosRev + nNegIrrev) speye(nNegRev)*-100 sparse(nNegRev, nPosIrrev + nPosRev*6 + nNegIrrev) speye(nNegRev) sparse(nNegRev,nEssRev*6)];    
%else
    nrRows1 = [sEye(negRevRxns,:) sparse(nNegRev,nYBlock + nPosIrrev + nPosRev*6 + nNegIrrev) speye(nNegRev)*-1 speye(nNegRev) sparse(nNegRev,nNegRev + nEssRev*6+nMetVars)];
    nrRows2 = [sparse(nNegRev,nRxns + nPosIrrev + nPosRev + nNegIrrev) speye(nNegRev)*-100 sparse(nNegRev,nPosIrrev + nPosRev*6 + nNegIrrev) speye(nNegRev) speye(nNegRev) speye(nNegRev) sparse(nNegRev,nEssRev*6+nMetVars)];
    nrRows = [nrRows1;nrRows2];
%end
erRows1 = [sEye(essRevRxns,:) sparse(nEssRev,nYBlock + nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev) speye(nEssRev)*-1 speye(nEssRev) sparse(nEssRev, nEssRev*4+nMetVars)];
erRows2 = [sparse(nEssRev,nRxns + nYBlock + nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev) speye(nEssRev) speye(nEssRev) speye(nEssRev)*-1 sparse(nEssRev, nEssRev*3+nMetVars)];
erRows3 = [sparse(nEssRev,nRxns + nYBlock + nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev) speye(nEssRev) sparse(nEssRev, nEssRev*2) speye(nEssRev)*-100 speye(nEssRev) sparse(nEssRev, nEssRev+nMetVars)];
erRows4 = [sparse(nEssRev,nRxns + nYBlock + nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev + nEssRev) speye(nEssRev) sparse(nEssRev, nEssRev) speye(nEssRev)*100 sparse(nEssRev, nEssRev) speye(nEssRev) sparse(nEssRev, nMetVars)];

%now the mets
if ~isempty(metData)
    %Order of rxn "on" vars: Ypi Ypr Yni Ynr. This is a bit messy, since they are not in the
    %same order as the reactions in the S matrix or metData. We therefore resort the vars in 
    %metData to match the order of the vars.
    srtMetData = sparse([metData(:,posIrrevRxns) metData(:,posRevRxns) metData(:,negIrrevRxns) metData(:,negRevRxns)]);
    
    %The mon vars come first followed by the mv1 vars
    metRows = [sparse(nMetabolMets,nRxns) -srtMetData sparse(nMetabolMets,  nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev + nEssRev*6) speye(nMetabolMets) speye(nMetabolMets)];
else
    metRows = [];
end

prob.a = [sRow;piRows;prRows1;prRows2;prRows3;prRows4;niRows;nrRows;erRows1;erRows2;erRows3;erRows4;metRows];
prob.A = prob.a;
prob.b = zeros(size(prob.A,1),1);
prob.c = [zeros(nRxns,1);rxnScores(posIrrevRxns)*-1;rxnScores(posRevRxns)*-1;rxnScores(negIrrevRxns)*-1;rxnScores(negRevRxns)*-1;zeros(nPosIrrev + nPosRev*6 + nNegIrrev + nNegRev*varsPerNegRev + nEssRev*6,1);ones(nMetabolMets,1)*-prodWeight;zeros(nMetabolMets,1)];
prob.lb = [milpModel.lb;zeros(nYBlock + nPosIrrev + 5*nPosRev,1);ones(nPosRev,1)*-100;zeros(nNegIrrev+nNegRev*varsPerNegRev+nEssRev*2,1);ones(nEssRev,1)*forceOnLim;zeros(nEssRev*2,1);ones(nEssRev,1)*-100; zeros(nMetVars,1)];
prob.lb(essIrrevRxns) = forceOnLim;%force flux for the ess irrev rxns - this is the only way these are handled
prob.ub = [milpModel.ub;ones(nYBlock,1);inf(nPosIrrev+nPosRev*2,1);ones(nPosRev,1);inf(3*nPosRev+nNegIrrev+varsPerNegRev*nNegRev+3*nEssRev,1);ones(nEssRev,1);inf(2*nEssRev,1);ones(nMetabolMets,1);inf(nMetabolMets,1)];

prob.vartype = [repmat('C', 1, nRxns + nPosIrrev + nPosRev), ...
                repmat('B', 1, nNegIrrev + nNegRev), ...
                repmat('C', 1, nPosIrrev + 2*nPosRev), ... 
                repmat('B', 1, nPosRev), ...
                repmat('C', 1, 3*nPosRev + nNegIrrev + varsPerNegRev*nNegRev + 3*nEssRev), ...
                repmat('B', 1, nEssRev), ...
                repmat('C', 1, 2*nEssRev + nMetVars)];

            
onoffVarInd = (1:(nPosIrrev + nPosRev + nNegIrrev + nNegRev)) + nRxns;
onoffPosIrrev = onoffVarInd(1:nPosIrrev);
onoffPosRev = onoffVarInd((1:nPosRev)+nPosIrrev);
onoffNegIrrev = onoffVarInd((1:nNegIrrev)+nPosIrrev+nPosRev);
onoffNegRev = onoffVarInd((1:nNegRev)+nPosIrrev+nPosRev+nNegIrrev);

metVarInd = (1:nMetabolMets) + (length(prob.vartype) - nMetVars);

if allowExcretion
    length(prob.b);
    length(milpModel.mets);
    prob.csense = [repmat('L', 1, length(milpModel.mets)), ...
                  repmat('E', 1, length(prob.b) - length(milpModel.mets))];
else
    prob.csense = '=';
end

params.intTol = 10^-8; %This value seems to work - making it smaller may slow the solver down.

prob.osense = 1; %minimize

res=optimizeProb(prob,params);


if ~checkSolution(res)
    if strcmp(res.origStat, 'TIME_LIMIT')
        EM='Time limit reached without finding a solution. Try increasing the TimeLimit parameter.';
    else
        EM='The problem is infeasible';
    end
    dispEM(EM);
end

%get the on/off vals
onoff = ones(nRxnsWithOnOff,1);%standard is on, i.e. all reactions not included in the problem + ess is on
onoff(posIrrevRxns) = res.full(onoffPosIrrev);
onoff(posRevRxns) = res.full(onoffPosRev);
onoff(negIrrevRxns) = res.full(onoffNegIrrev);
onoff(negRevRxns) = res.full(onoffNegRev);

onoff2 = zeros(nRxnsWithOnOff,1);%standard is off, i.e. all reactions not included in the problem + ess is off
onoff2(posIrrevRxns) = res.full(onoffPosIrrev);
onoff2(posRevRxns) = res.full(onoffPosRev);
onoff2(negIrrevRxns) = res.full(onoffNegIrrev);
onoff2(negRevRxns) = res.full(onoffNegRev);


%investigate a bit
%problematic = onoff < 0.99 & onoff > 0.01;
%onoff(problematic)

%so, it is only positive that are problematic. Likely because the 0.1 is too large. Test to change to 0.01
%unique(onoff(onoff < 0.99 & onoff > 0.01))

%Get all reactions used in the irreversible model
%The reaction score check is for metabolomics - we add those reactions as 
%positive even though the score is 0. And, we don't want them included in the 
%results.
deletedRxns=(onoff < 0.5).' & (rxnScores ~= 0).'; 
turnedOnRxns=(onoff2 >= 0.5).' & (rxnScores ~= 0).';

fluxes = res.full(1:nRxns);

%extract the met data
metProduction = logical(round(res.full(metVarInd)));

end
