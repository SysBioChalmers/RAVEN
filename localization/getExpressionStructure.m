function experiment=getExpressionStructure(fileName)
% getExpressionStructure
%   Loads a representation of an experiment from an Excel file (see
%   comments further down)
%
%   fileName            an Excel representation on an experiment
%
%   experiment          an experiment structure
%       data            matrix with expression values
%       orfs            the corresponding ORFs
%       experiments     the titles of the experiments
%       boundNames      reaction names for the bounds
%       upperBoundaries matrix with the upper bound values
%       fitNames        reaction names for the measured fluxes
%       fitTo           matrix with the measured fluxes
%
%   A very common data set when working with genome-scale metabolic models
%   is that you have measured fermentation data, gene expression data,
%   and some different 'bounds' (for example different carbon sources
%   or genes that are knocked out) in a number of conditions. This function
%   reads an Excel representation of such an experiment.
%   The Excel file must contain three sheets, 'EXPRESSION', 'BOUNDS',
%   'FITTING'. Below are some examples to show how they should be
%   formatted:
%
%   -EXPRESSION
%       ORF	dsm_paa	wisc_paa
%       Pc00e00030	79.80942723	78.14755338
%   Shows the expression of the gene Pc00e00030 under two different
%   conditions (in this case a DSM strain and a Wisconsin strain of P.
%   chrysogenum with PSS in the media)
%
%   -BOUNDS
%       Fixed Upper	dsm_paa	wisc_paa
%       paaIN	0.1	0.2
%   The upper bound for the reaction paaIN should be 0.1 for the first
%   condition and 0.2 for the second
%
%   -FITTING
%       Fit to	dsm_paa	wisc_paa
%       co2OUT	2.85	3.05
%       glcIN   1.2     0.9
%   The measured fluxes for CO2 production and glucose uptake for the two
%   conditions. The model(s) can later be fitted to match these values as
%   good as possible.
%
% Usage: experiment=getExpressionStructure(fileName)

[type, sheets]=xlsfinfo(fileName);

%Check if the file is a Microsoft Excel Spreadsheet
if ~strcmp(type,'Microsoft Excel Spreadsheet')
    EM='The file is not a Microsoft Excel Spreadsheet';
    dispEM(EM);
end

%Check that all sheets are present and saves the index of each
exprSheet=find(strcmp('EXPRESSION', sheets));
boundSheet=find(strcmp('BOUNDS', sheets));
fitSheet=find(strcmp('FITTING', sheets));

if length(exprSheet)~=1 || length(boundSheet)~=1 || length(fitSheet)~=1
    EM='Not all required spreadsheets are present in the file';
    dispEM(EM);
end

%Load the expression data
[discard,dataSheet]=xlsread(fileName,exprSheet);
experiment.data=discard;
experiment.orfs=dataSheet(2:size(dataSheet,1),1);
experiment.experiments=dataSheet(1,2:size(dataSheet,2));

%Loads the maximal boundaries
[discard,dataSheet]=xlsread(fileName,boundSheet);
experiment.boundNames=dataSheet(2:size(dataSheet,1),1);
experiment.upperBoundaries=discard;

%Loads the experimental data to fit to
[discard,dataSheet]=xlsread(fileName,fitSheet);
experiment.fitNames=dataSheet(2:size(dataSheet,1),1);
experiment.fitTo=discard;

%Check to see that the dimensions are correct
if length(experiment.orfs)~=size(experiment.data,1) || (length(experiment.experiments)~=size(experiment.data,2) && ~isempty(experiment.data))
    EM='The expression data does not seem to be formated in the expected manner';
    dispEM(EM);
end
end
