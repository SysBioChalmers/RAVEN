function geneScoreStructure=parseScores(inputFile,predictor)
% parseScores
%   Parses the output from a predictor to generate the geneScoreStructure.
%
%   inputFile	a file with the output from the predictor
%   predictor   the predictor that was used. 'tsv' for tab-separated values
%               where the name of the compartments in the first row and each
%               row after that correspond to a gene. 'wolf' for 
%               WoLFPSORT. (opt, default 'tsv')
%
%   The function normalizes the scores so that the best score for each gene
%   is 1.0.
%
%   geneScoreStructure  a structure to be used in predictLocalization
%
%   Usage: geneScoreStructure=parseScores(inputFile,predictor,normalize)
%
%   Rasmus Agren, 2013-08-01

if nargin<2
    predictor='tsv';
end

fid=fopen(inputFile,'r');

if fid<1
   dispEM('Could not open file');  
end

if strcmpi(predictor,'wolf')
   A=textscan(fid,'%s','Delimiter','\n','CommentStyle','#'); 
   
   %Each element should be for one gene, but some of them are on the form
   %"Pc20g11350: treating 9 X's as Glycines". Those should be removed.
   I=~cellfun(@any,strfind(A{1},'treating'));
   
   B=regexp(A{1}(I),' ','split');
   
   %Reserve space for stuff
   geneScoreStructure.compartments={};
   geneScoreStructure.scores=[]; %Don't know number of comps yet
   geneScoreStructure.genes=cell(numel(B),1);
   
   %Parsing is a bit cumbersome as ', ' is used as a delimiter in some cases
   %and ' ' in others. Use strrep to get rid of ','.
   for i=1:numel(B)
        b=strrep(B{i},',','');
        geneScoreStructure.genes{i}=b{1};
        
        %Then go through the compartments and add new ones as they are
        %found
        for j=2:2:numel(b)-1
        	[crap J]=ismember(b(j),geneScoreStructure.compartments);
            
            %Add new compartment if it doesn't exist
            if J==0
               geneScoreStructure.compartments=[geneScoreStructure.compartments;b(j)];
               J=numel(geneScoreStructure.compartments);
               geneScoreStructure.scores=[geneScoreStructure.scores zeros(numel(B),1)];
            end
            
            geneScoreStructure.scores(i,J)=str2double(b(j+1));
        end
   end
end

%Check if there are duplicate genes
[crap J K]=unique(geneScoreStructure.genes);

if numel(J)~=numel(K)
   dispEM('There are duplicate genes in the input file',false);
   geneScoreStructure.genes=geneScoreStructure.genes(J);
   geneScoreStructure.scores=geneScoreStructure.scores(J,:);   
end

%Normalize
I=max(geneScoreStructure.scores,[],2);
geneScoreStructure.scores=bsxfun(@times, geneScoreStructure.scores, 1./I);

fclose(fid);