function FSEOF(model,biomassRxn,targetRxn,iterations,coefficient,outputFile)
% FSEOF: implements the algorithm of Flux Scanning based on Enforced Objective Flux.
%
%   model           a model structure
%   biomassRxn      string with reaction ID of the biomass formation or
%                   growth reaction
%   targetRxn       string with reaction ID of target reaction
%   iterations      double indicating number of iterations (opt, default 10)
%   coefficient     double indicating ratio of optimal target reaction
%                   flux, must be less than 1 (opt, default 0.9)
%   outputFile      string with output filename (opt, default prints to
%                   command window)
%
%OUTPUTS
%   This function output an tab-delimited file or prints to command window.
%
%   Usage: FSEOF(model,biomassRxn,targetRxn,iterations,coefficient,outputFile)
%
%   Eduard Kerkhoven, 2017-06-02

output = 0;
if nargin<4
  iterations = 10;
  coefficient = 0.9;
end

if nargin <5
	coefficient = 0.9;
end

if nargin == 6
	output = 1;
end

%Find out the maximum theoretical yield of target reaction
model = setParam(model,'obj',targetRxn,1);
sol=solveLP(model,1);
targetMax = abs(sol.f*coefficient);   % 90 percent of the theoretical yield

model = setParam(model,'obj',biomassRxn,1);

fseof.results=zeros(length(model.rxns),iterations);
fseof.target=zeros(length(model.rxns),1);
rxnDirection=zeros(length(model.rxns),1);

%Enforce objective flux iteratively
for i=1:iterations
	n = i*targetMax/iterations;
	model = setParam(model,'lb',targetRxn,n);
  
  sol=solveLP(model,1);
  
	fseof.results(:,i)=sol.x;

	%Loop through all fluxes and identify the ones that increased upon the enforced objective flux
	for j=1:length(fseof.results)
  	if fseof.results(j,1) > 0   %Check the positive fluxes
  	
  		if i == 1   %The initial round
  			rxnDirection(j,1) = 1;
  			fseof.target(j,1) = 1;
			else
				
				if (fseof.results(j,i) > fseof.results(j,i-1)) & fseof.target(j,1)
  				fseof.target(j,1) = 1;
  			else
  				fseof.target(j,1) = 0;
  			end
  		end

  	elseif fseof.results(j,1) < 0 %Check the negative fluxes
  	
  		if i == 1   %The initial round
  			rxnDirection(j,1) = -1;
  			fseof.target(j,1) = 1;
			else
				if (fseof.results(j,i) < fseof.results(j,i-1)) & fseof.target(j,1)
  				fseof.target(j,1) = 1;
  			else
  				fseof.target(j,1) = 0;
  			end
  		end
  	
		end

	end
end

%Generating output
formatSpec = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n';
if output == 1    %Output to a file
	fid = fopen(outputFile,'w');
	fprintf(fid,formatSpec,'Slope','rowID','Enzyme ID','Enzyme Name','Subsystem','Direction','Gr Rule');
else              %Output to screen
	fprintf(formatSpec,'Slope','rowID','Enzyme ID','Enzyme Name','Subsystem','Direction','Gr Rule');
end

for num=1:length(fseof.target)
	if fseof.target(num,1) == 1
		A0 = num2str(abs(fseof.results(num,iterations)-fseof.results(num,1))/abs(targetMax-targetMax/iterations)); %Slope calculation
		A1 = num2str(num);                                  %row ID
		A2 = char(model.rxns(num));                         %enzyme ID
		A3 = char(model.rxnNames(num));                     %enzyme Name
		A4 = char(model.subSystems(num));                   %Subsystem
		A5 = num2str(model.rev(num)*rxnDirection(num,1));   %reaction Dirction
		A6 = char(model.grRules(num));                      %Gr Rule
		if output == 1    %Output to a file
			fprintf(fid,formatSpec,A0,A1,A2,A3,A4,A5,A6);
		else              %Output screen
			fprintf(formatSpec,A0,A1,A2,A3,A4,A5,A6);
		end
	end
end

if output == 1    %Output to a file
	fclose(fid);
end
end
