function [S, mets, badRxns, reversible]=constructS(equations,mets,rxns)
% constructS
%   Constructs a stoichiometric matrix from a cell array of equations
%
%   equations   cell array of equations on the form 'A + 2 B <=> or => 3 C'
%   mets        cell array of metabolites. All metabolites in the equations
%               must be present in the list (opt, default generated from the equations)
%   rxns        cell array of reaction ids. This is only used for printing
%               reaction ids instead of equations in warnings/errors (opt,
%               default [])
%
%
%   S           the resulting stoichiometric matrix
%   mets        cell array with metabolites that corresponds to the order in
%               the S matrix
%   badRxns     boolean vector with the reactions that have one or more
%               metabolites as both substrate and product. An example would be
%               the phosphotransferase ATP + ADP <=> ADP + ATP. In the
%               stoichiometric matrix this equals to an empty reaction which
%               can be problematic
%   reversible  boolean vector with true if the equation is reversible
%
%   Usage: [S mets badRxns reversible]=constructS(equations,mets)
%
%   Rasmus Agren, 2014-01-08
%

badRxns=false(numel(equations),1);

%Check that no equations are too short to have reversibility data
I=cellfun(@numel,equations);
I=find(I<4,1);
if any(I)
    if isempty(rxns)
        EM=['The following equation does not have reversibility data: ' equations{I} ];
        dispEM(EM);
    else
        EM=['The reaction ' rxns{I} ' does not have reversibility data'];
        dispEM(EM);
    end
end

%Makes life a little easier
equations=strtrim(equations);
equations=fixEquations(equations);

if nargin<2
    mets=parseRxnEqu(equations);
end
if nargin<3
    rxns=[];
end

%Get which reactions are reversible
reversible=cellfun(@any,strfind(equations,' <=> '));

%Make them all reversible. This is not all that neat, but nevermind
equations=strrep(equations,' => ',' <=> ');

%Replace the the plus signs with some weird character that will be used for
%parsing
equations=strrep(equations,' + ', '?');

%Generate the stoichiometric matrix
S=zeros(numel(mets),numel(equations));

%Loop through the equations and add the info to the S matrix
for i=1:numel(equations)
    %Start by finding the position of the (=> or <=>)
    arrowIndex=strfind(equations{i},' <=> ');

    if numel(arrowIndex)~=1
        if isempty(rxns)
            EM=['The following equation does not have reversibility data: ' equations{i} ];
            dispEM(EM);
        else
            EM=['The reaction ' rxns{i} ' does not have reversibility data'];
            dispEM(EM);
        end
    end

    reactants=regexp(equations{i}(1:arrowIndex-1),'?','split');
    products=regexp(equations{i}(arrowIndex+5:end),'?','split');

    %If the splitting character is at the end (if exchange rxns), then an
    %empty string will exist together with the real ones. Remove it
    reactants(cellfun(@isempty,reactants))=[];
    products(cellfun(@isempty,products))=[];

    %A vector where an element is -1 is the corresponding metabolite is a
    %reactant and 1 if it's a product
    multiplyWith=[ones(numel(reactants),1)*-1; ones(numel(products),1)];

    metabolites=[reactants products];

    %Now loop through the reactants and see if the metabolite has a coefficient
    %(it will look as 'number name')
    for j=1:numel(metabolites)
        space=strfind(metabolites{j},' ');

        if isempty(space)
            %No coefficient
            coeff=1;
            name=metabolites{j};
        else
            coeff=str2double(metabolites{j}(1:space(1)));

            %If it was not a coefficiant
            if isnan(coeff)
               coeff=1;
               name=metabolites{j};
            else
               name=metabolites{j}(space+1:end);
            end
        end

        %Find the name in the mets list
        %[a b]=ismember(name,mets);
        b=find(strcmp(name,mets),1);

        if any(b)
            %Check if the metabolite already participates in this reaction
            if S(b,i)~=0
               badRxns(i)=true;
            end
            S(b,i)=S(b,i)+coeff*multiplyWith(j);
        else
            if isempty(rxns)
                EM=['Could not find metabolite ' name ' in metabolite list'];
                dispEM(EM);
            else
                EM=['The metabolite "' name '" in reaction ' rxns{i} ' was not found in the metabolite list'];
                dispEM(EM);
            end
        end
    end
end
S=sparse(S);
end

function equ=fixEquations(equ)
%If the equation starts with "=>" or "<=>" then add a space again. This is
%an alternative way to represent uptake reactions. The opposite way for
%producing reactions
equ=equ(:);
for i=1:numel(equ)
   if strcmp(equ{i}(1:2),'=>') || strcmp(equ{i}(1:3),'<=>')
       equ{i}=[' ' equ{i}];
   else
       if strcmp(equ{i}(end-1:end),'=>') || strcmp(equ{i}(end-2:end),'<=>')
            equ{i}=[equ{i} ' '];
       end
   end
end
end
