function equationString=buildEquation(mets,stoichCoeffs,isrev)
% buildEquation
%   Construct single equation string for a given reaction
%
%   mets            cell array with metabolites involved in the reaction.
%   stoichCoeffs    vector with corresponding stoichiometric coeffs.
%   isrev           logical indicating if the reaction is or not
%                   reversible.
%
%   equationString  equation as a string
%
%    Usage: equationString=buildEquation(mets,stoichCoeffs,isrev)

mets=convertCharArray(mets);
if ~isnumeric(stoichCoeffs)
    EM = 'stoichCoeffs must be a numeric vector';
    dispEM(EM);
elseif ~islogical(isrev)
    EM = 'isrev must be a logical';
    dispEM(EM);
elseif length(mets) ~= length(stoichCoeffs)
    EM = 'lengths of mets and stoichCoeffs should be the same';
    dispEM(EM);
end

%Reactant half:
reactants       = mets(stoichCoeffs<0);
reactantsCoeffs = abs(stoichCoeffs(stoichCoeffs<0));
reactantsEqn    = concatenateEquation(reactants,reactantsCoeffs);

%Product half:
products       = mets(stoichCoeffs>0);
productsCoeffs = abs(stoichCoeffs(stoichCoeffs>0));
productsEqn    = concatenateEquation(products,productsCoeffs);

%Full equation:
if isrev
    equationString = [reactantsEqn ' <=> ' productsEqn];
else
    equationString = [reactantsEqn ' => ' productsEqn];
end

end

function eqn = concatenateEquation(mets,stoichCoeffs)
%This function concatenates metabolites and stoich. coefficients to form
%either the left or right side of the rxn equation
eqn = '';
for i = 1:length(stoichCoeffs)
    if i == 1
        plusString='';
    else
        plusString=' + ';
    end
    stoich = stoichCoeffs(i);
    if stoich == 1
        stoich = '';
    else
        stoich = [num2str(stoich) ' '];
    end
    eqn = [eqn plusString stoich mets{i}];
end
end
