function [elements, useMat, exitFlag, MW]=parseFormulas(formulas, noPolymers,isInchi,ignoreRX)
% parseFormulas
%   Gets the elemental composition from formulas
%
%   formulas      a cell array with formulas
%   noPolymers    assume that all polymers consist of one element.
%                 Corresponds to counting everything between (...)n as
%                 n being equal to one. Only one set of parentheses
%                 is allowed. If this is false then polymers are returned as
%                 "Could not parse formula" (opt, default false)
%   isInchi       true if the formulas are in the InChI format (opt,
%                 default false)
%   ignoreRX      ignore R-groups and bound protein. This can be useful since they
%                 are often used only as intermediates (opt, default false)
%
%   elements
%       abbrevs   cell array with abbreviations for all used elements
%       names     cell array with the names for all used elements
%   useMat        MxN matrix with the number of atoms for each formula (M) and each
%                 element (N)
%   exitFlag      array with the exit flags:
%                 1=  Sucessful parsing
%                 0=  No formula found
%                 -1= Could not parse formula
%   MW            predicted molecular weight (g/mol). This is only returned
%                 for formulas which can be sucessfully parsed, and its
%                 calculation doesn't affect the exitFlag variable. NaN is
%                 returned if the weight couldn't be calculated
%   
%   Usage: [elements, useMat, exitFlag, MW]=
%               parseFormulas(formulas, noPolymers,isInchi,ignoreRX)

if nargin<2
    noPolymers=false;
end
if nargin<3
    isInchi=false;
end
if nargin<4
    ignoreRX=false;
end

elements.abbrevs={'C', 'N', 'O', 'S', 'P', 'H', 'He', 'Li', 'Be', 'B', 'F', 'Ne', 'Na', 'Mg', 'Al',...
    'Si', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',...
    'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc',...
    'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce',...
    'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',...
    'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',...
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',...
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'R', 'X'}';
elements.names={'carbon', 'nitrogen', 'oxygen', 'sulfur', 'phosphorus', 'hydrogen', 'helium', 'lithium', 'beryllium', 'boron',...
    'fluorine', 'neon', 'sodium', 'magnesium', 'aluminum,', 'silicon',...
    'chlorine', 'argon', 'potassium', 'calcium', 'scandium', 'titanium', 'vanadium',...
    'chromium', 'manganese', 'iron', 'cobalt', 'nickel', 'copper', 'zinc', 'gallium', 'germanium',...
    'arsenic', 'selenium', 'bromine', 'krypton', 'rubidium', 'strontium', 'yttrium', 'zirconium',...
    'niobium', 'molybdenum', 'technetium', 'ruthenium', 'rhodium', 'palladium', 'silver', 'cadmium',...
    'indium', 'tin', 'antimony', 'tellurium', 'iodine', 'xenon', 'cesium', 'barium', 'lanthanum',...
    'cerium', 'praseodymium', 'neodymium', 'promethium', 'samarium', 'europium', 'gadolinium',...
    'terbium', 'dysprosium', 'holmium', 'erbium', 'thulium', 'ytterbium', 'lutetium', 'hafnium',...
    'tantalum', 'tungsten', 'rhenium', 'osmium', 'iridium', 'platinum', 'gold', 'mercury',...
    'thallium', 'lead', 'bismuth', 'polonium', 'astatine', 'radon', 'francium', 'radium',...
    'actinium', 'thorium', 'protactinium', 'uranium', 'neptunium', 'plutonium', 'americium',...
    'curium', 'berkelium', 'californium', 'einsteinium', 'fermium', 'mendelevium', 'nobelium',...
    'lawrencium', 'rutherfordium', 'dubnium', 'seaborgium', 'bohrium', 'hassium', 'meitnerium',...
    'darmstadtium', 'roentgenium', 'copernicium', 'generic group', 'bound protein'}';

EWs=[12.0107 14.0067 15.9994 32.065 30.973762 1.00794 4.002602 6.941 9.012182 10.811 18.9984032 ...
    20.1797 22.98976928 24.305 26.9815386 28.0855 35.453 39.948 39.0983 40.078 44.955912 47.867 50.9415 ...
    51.9961 54.938045 55.845 58.933195 58.6934 63.546 65.39 69.723 72.64 74.9216 78.96 79.904 83.798 ...
    85.4678 87.62 88.90585 91.224 92.906 95.94 97.9072 101.07 102.905 106.42 107.8682 112.411 114.818 ...
    118.71 121.76 127.6 126.904 131.293 132.9054519 137.327 138.90547 140.116 140.90765 144.242 144.9127 ...
    150.36 151.964 157.25 158.92535 162.5 164.93 167.259 168.93421 173.04 174.967 178.49 180.94788 183.84 ...
    186.207 190.23 192.217 195.084 196.966569 200.59 204.3833 207.2 208.9804 208.9824 209.9871 222.0176 ...
    223.0197 226.0254 227.0277 232.03806 231.03588 238.02891 237.0482 244.0642 243.0614 247.0704 247.0703 ...
    251.0796 252.083 257.0951 258.0984 259.101 262.1097 261.1088 262 266 264 277 268 271 272 nan nan nan]';

%Set the EWs of these groups to 0
if ignoreRX==true
    EWs(end-1:end)=0;
end

useMat=zeros(numel(formulas),numel(elements.abbrevs));

exitFlag=zeros(numel(formulas),1);

%"H+", if parsed from InChI code, would have the composition "p+1". Change
%to fit with how the rest of the compounds are written
formulas=strrep(formulas,'p+1','H+');

%Ignore charge to make the parsing easier
formulas=strrep(formulas,'+','');
formulas=strrep(formulas,'-','');

%Loop through each formula
for i=1:numel(formulas)
    if ~isempty(formulas{i})
        sucess=false; %To see if it works
        formula=formulas{i};
        
        %If it's an InChI code. The composition is found between the first
        %and the second "/" For some simple molecules such as salts only
        %the first "/" is present
        if isInchi==true
            S=regexp(formula,'/','split');
            if numel(S)>=2
                formula=S{2};
            else
                formula='';
            end
        end
        %Only look at what's between the parantheses (polymers are not
        %supported in InChI)
        if isInchi==false
            LP=strfind(formula,'(');
            RP=strfind(formula,')n');
            %            if numel(LP)>1 || numel(RP)>1
            %               exitFlag(i)=-1; continue;
            %            end if numel(LP)>1 || numel(RP)>1
            %               exitFlag(i)=-1; continue;
            %            end
            if numel(LP)==1 && numel(RP)==1
                %This means that the polymer should be regarded as a
                %monomer
                if noPolymers==true
                    %This means that there are one set of parantheses
                    formula=strrep(formula,'(','');
                    formula=strrep(formula,')n','');
                else
                    %This means that polymers should be ignored
                    exitFlag(i)=-1;
                    continue;
                end
            else
                if ~isempty(LP) || ~isempty(RP)
                    exitFlag(i)=-1;
                    continue;
                end
            end
        end
        
        %Get the indexes of the numeric (or ".") characters
        nonNumeric=false(numel(formula),1);
        nonNumeric(regexp(formula,'[^0-9.]'))=true;
        
        %Get the indexes of the upper characters (since each element starts
        %with an upper character)
        upperI=isstrprop(formula,'upper');
        upperX=find(upperI);
        
        for j=1:numel(upperX)
            %The first case is when it's the last character. Then the
            %coefficient must be 1
            isLast=false;
            if upperX(j)==numel(formula)
                coeff=1;
                element=formula(upperX(j));
                isLast=true;
            end
            
            if isLast==false
                %The second case is when the following character is a
                %character
                if nonNumeric(upperX(j)+1)
                    %Is it a new element?
                    if upperI(upperX(j)+1)
                        %New element, that means that the coefficient was 1
                        %and that the element was only one character
                        coeff=1;
                        element=formula(upperX(j));
                    else
                        %This means that it's an element with two
                        %characters
                        if j==numel(upperX)
                            if upperX(j)<numel(formula)-1
                                coeff=str2double(formula(upperX(j)+2:end));
                            else
                                coeff=1;
                            end
                        else
                            %Check if there is a number or a new element
                            %after it
                            if nonNumeric(upperX(j)+2)==true
                                coeff=1;
                            else
                                coeff=str2double(formula(upperX(j)+2:upperX(j+1)-1));
                            end
                        end
                        element=formula(upperX(j):upperX(j)+1);
                    end
                else
                    %Then it is a numeral
                    if j==numel(upperX)
                        coeff=str2double(formula(upperX(j)+1:end));
                    else
                        coeff=str2double(formula(upperX(j)+1:upperX(j+1)-1));
                    end
                    element=formula(upperX(j));
                end
            end
            
            %Find the element
            I=strcmp(element,elements.abbrevs);
            if any(I)
                if ~isnan(coeff)
                    useMat(i,I)=useMat(i,I)+coeff;
                    sucess=true;
                else
                    break;
                end
            else
                break;
            end
        end
        if sucess==false
            useMat(i,:)=0; %Reset for this formula
            exitFlag(i)=-1;
        else
            exitFlag(i)=1;
        end
    end
end

%Remove the elements which are never used
I=~any(useMat);
useMat(:,I)=[];
elements.abbrevs(I)=[];
elements.names(I)=[];
EWs(I)=[];

%Calcluate the weight contribution of each element in each formula. Note
%that this will give NaNs for all formulas if R or X groups are in EWs,
%since 0*NaN is still NaN. Therefore only use elements with known mass
if nargout>3
    P=bsxfun(@times,useMat(:,~isnan(EWs)),EWs(~isnan(EWs)).');
    MW=sum(P,2);
    
    %Then remove the calculations for elements with unknown mass
    I=find(useMat(:,isnan(EWs)));
    MW(I)=nan;
    MW(exitFlag~=1)=nan;
end
end
