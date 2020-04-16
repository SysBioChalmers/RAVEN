function GSS = parseScores(inputFile, predictor)
% parseScores
%	Parse the output from a predictor to generate the GSS
%
%	Input:
%	inputFile	a file with the output from the predictor
%	predictor	the predictor that was used. 'wolf' for WoLF PSORT, 'cello'
%               for CELLO, 'deeploc' for DeepLoc (opt, default 'wolf')
%
%	Output:
%	GSS         a gene scoring structure to be used in predictLocalization
%
%	The function normalizes the scores so that the best score for each gene
%	is 1.0.
%
%	Usage: GSS = parseScores(inputFile, predictor)
%
%	Simonas Marcisauskas, 2019-11-13
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
    %"Pc20g11350: treating 9 X's as Glycines". Those should be removed
    I=~cellfun(@any,strfind(A{1},'treating'));
    
    B=regexp(A{1}(I),' ','split');
    
    %Reserve space for stuff
    GSS.compartments={};
    GSS.scores=[]; %Do not know number of comps yet
    GSS.genes=cell(numel(B),1);
    
    %Parsing is a bit cumbersome as ', ' is used as a delimiter in some
    %cases and ' ' in others. Use strrep to get rid of ','
    for i=1:numel(B)
        b=strrep(B{i},',','');
        GSS.genes{i}=b{1};
        
        %Then go through the compartments and add new ones as they are
        %found
        for j=2:2:numel(b)-1
            [~, J]=ismember(b(j),GSS.compartments);
            
            %Add new compartment if it does not exist
            if J==0
                GSS.compartments=[GSS.compartments;b(j)];
                J=numel(GSS.compartments);
                GSS.scores=[GSS.scores zeros(numel(B),1)];
            end
            
            GSS.scores(i,J)=str2double(b(j+1));
        end
    end
elseif strcmpi(predictor,'cello')
    fid=fopen(inputFile,'r');
    %Read the title line and fetch the list of compartments
    tline = fgetl(fid);
    tline=regexprep(tline,'^.+#Combined:\t','');
    tline=regexprep(tline,'\t#Most-likely-Location.+','');
    GSS.compartments=transpose(regexp(tline,'\t','split'));
    
    %Now iterate through the following lines in the file. Each row
    %corresponds to one gene and it consists of the scores for
    %compartments. Gene name is in the end of each line
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
        GSS.scores(row,:)=str2double(tline(1:numel(GSS.compartments)));
        GSS.genes{row,1}=tline{1,end};
    end
elseif strcmpi(predictor,'deeploc')
    fid=fopen(inputFile,'r');
    %Read the title line and fetch the list of compartments
    tline = fgetl(fid);
    GSS.compartments=regexp(tline,'\t','split');
    GSS.compartments=GSS.compartments(3:end);
    
    %Now iterate through the following lines in the file. Each row
    %corresponds to one gene and it consists of the scores for
    %compartments. Gene name is in the end of each line
    row=0;
    while 1
        row=row+1;
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        tline=regexp(tline,'\t','split');
        GSS.scores(row,:)=str2double(tline(3:end));
        GSS.genes{row,1}=tline{1,1};
    end
end

%Check if there are duplicate genes
[~, J, K]=unique(GSS.genes);

if numel(J)~=numel(K)
    EM='There are duplicate genes in the input file';
    dispEM(EM,false);
    GSS.genes=GSS.genes(J);
    GSS.scores=GSS.scores(J,:);
end

%Normalize
I=max(GSS.scores,[],2);
GSS.scores=bsxfun(@times, GSS.scores, 1./I);

fclose(fid);
end
