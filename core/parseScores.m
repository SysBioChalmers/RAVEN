function geneScoreStructure=parseScores(inputFile,predictor)
% parseScores
%   Parses the output from a predictor to generate the geneScoreStructure.
%
%   inputFile	a file with the output from the predictor
%   predictor   the predictor that was used. 'wolf' for WoLF PSORT, 'cello'
%               for CELLO. (opt, default 'wolf')
%
%   The function normalizes the scores so that the best score for each gene
%   is 1.0.
%
%   geneScoreStructure  a structure to be used in predictLocalization
%
%   Usage: geneScoreStructure=parseScores(inputFile,predictor,normalize)
%
%   Simonas Marcisauskas, 2016-11-15 - added compatibility for CELLO v2.5
%

if nargin<2
    predictor='wolf';
end

fid=fopen(inputFile,'r');

if fid<1
   EM='Could not open file';
   dispEM(EM);
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
        	[~, J]=ismember(b(j),geneScoreStructure.compartments);

            %Add new compartment if it doesn't exist
            if J==0
               geneScoreStructure.compartments=[geneScoreStructure.compartments;b(j)];
               J=numel(geneScoreStructure.compartments);
               geneScoreStructure.scores=[geneScoreStructure.scores zeros(numel(B),1)];
            end

            geneScoreStructure.scores(i,J)=str2double(b(j+1));
        end
   end
else if strcmpi(predictor,'cello')
    fid=fopen(inputFile,'r');
    % Reading the title line and fetching the list of compartments;
    tline = fgetl(fid);
    tline=regexprep(tline,'^.+#Combined:\t','');
    tline=regexprep(tline,'\t#Most-likely-Location.+','');
    geneScoreStructure.compartments=transpose(regexp(tline,'\t','split'));

    % Now iterating through the following lines in the file. Each row
    % corresponds to one gene and it consists of the scores for
    % compartments. Gene name is in the end of each line;
    row=0;
    while 1
        row=row+1;
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        tline=regexprep(tline,'^.+:\t','');
        tline=regexprep(tline,' .+','');
        tline=regexp(tline,'\t','split');
        geneScoreStructure.scores(row,:)=str2double(tline(1:numel(geneScoreStructure.compartments)));
        geneScoreStructure.genes{row,1}=tline{1,end};
        end
    end
end

%Check if there are duplicate genes
[crap J K]=unique(geneScoreStructure.genes);

if numel(J)~=numel(K)
   EM='There are duplicate genes in the input file';
   dispEM(EM,false);
   geneScoreStructure.genes=geneScoreStructure.genes(J);
   geneScoreStructure.scores=geneScoreStructure.scores(J,:);
end

%Normalize
I=max(geneScoreStructure.scores,[],2);
geneScoreStructure.scores=bsxfun(@times, geneScoreStructure.scores, 1./I);

fclose(fid);
