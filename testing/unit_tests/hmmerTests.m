%run this test case with the command
%results = runtests('hmmerTests.m')
function tests = hmmerTests
tests = functiontests(localfunctions);
end

function testHmmer(testCase)
%Get the directory for RAVEN Toolbox
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));

if isunix
    if ismac
        binEnd='.mac';
    else
        binEnd='';
    end
elseif ispc
    binEnd='.exe';
else
    dispEM('Unknown OS, exiting.')
    return
end

%Create empty structures needed for HMMER results
actHmmResult.genes = {};
actHmmResult.scores = [];

%Create structures that contain expected HMMER results
expHmmResult.genes = {'sp|P41947|MEL6_YEASX','sp|P41946|MEL5_YEASX', 'sp|P41945|MEL2_YEASX', 'sp|P04824|MEL1_YEASX'};
expHmmResult.scores = [0, 0, 0, 0];

%Generate temporary names for working directory and outFile
tmpDIR=tempname;
outFile=tempname;

%Run HMMER multi-threaded to use all logical cores assigned to MATLAB
cores = evalc('feature(''numcores'')');
cores = strsplit(cores, 'MATLAB was assigned: ');
cores = regexp(cores{2},'^\d*','match');
cores = cores{1};

%Create a temporary folder and copy multi-FASTA file there
[~, ~]=system(['mkdir "' tmpDIR '"']);

sourceDir = fileparts(which(mfilename));
copyfile(fullfile(sourceDir,'test_data','yeast_galactosidases.fa'),tmpDIR);

%Train a hidden Markov model
[~, ~]=system(['"' fullfile(ravenPath,'software','hmmer',['hmmbuild' binEnd]) '" --cpu "' num2str(cores) '" "' fullfile(tmpDIR,'yeast_galactosidases.hmm') '" "' fullfile(tmpDIR,'yeast_galactosidases.fa') '"']);

%Run a homology search against the newly-trained HMM
[~, output]=system(['"' fullfile(ravenPath,'software','hmmer',['hmmsearch' binEnd]) '" --cpu "' num2str(cores) '" "' fullfile(tmpDIR,'yeast_galactosidases.hmm') '" "' fullfile(tmpDIR,'yeast_galactosidases.fa') '"']);

%Save the output to a file
fid=fopen(outFile,'w');
fwrite(fid,output);
fclose(fid);


%Parse the results
geneCounter=1;

fid=fopen(outFile,'r');
beginMatches=false;
while 1
    %Get the next line
    tline = fgetl(fid);
    
    %Abort at end of file
    if ~ischar(tline)
        break;
    end
    
    if and(beginMatches,strcmp(tline,'  ------ inclusion threshold ------'))
        break;
    end
    
    if beginMatches==false
        %This is how the listing of matches begins
        if any(strfind(tline,'E-value '))
            %Read one more line that is only padding
            tline = fgetl(fid);
            beginMatches=true;
        end
    else
        %If matches should be read
        if ~strcmp(tline,'   [No hits detected that satisfy reporting thresholds]') && ~isempty(tline)
            elements=regexp(tline,' ','split');
            elements=elements(cellfun(@any,elements));
            
            actHmmResult.scores(geneCounter)=str2double(elements{1});
            actHmmResult.genes{geneCounter}=elements{9};
            
            geneCounter=geneCounter+1;
        else
            break;
        end
    end
end
fclose(fid);


%Remove the old tempfiles
delete([outFile '*']);

%Remove temporary folder, since testing is finished
[~, ~]=system(['rm "' tmpDIR '" -r']);


%Check if HMMER results match the expected ones
verifyEqual(testCase,actHmmResult,expHmmResult);

%Modify actual HMMER results structure and check if test fails
actHmmResult.score(2)=1;
verifyNotEqual(testCase,actHmmResult,expHmmResult);
end
