function model=getRxnsFromMetaCyc(metacycPath,keepTransportRxns,keepUnbalanced,keepUndetermined)
% getRxnsFromMetaCyc
%   Retrieves reactions information from MetaCyc database
%
%   Input:
%   metacycPath         if metaCycRxns.mat is not in the RAVEN\external\metacyc
%                       directory, this function will attempt to build it by
%                       reading info from a local dump of MetaCyc database,
%                       metacycPath is the path to the MetaCyc data files
%   keepTransportRxns   include transportation reactions, which often have identical
%                       reactants and products that turn to be all-zero columns in
%                       the S matrix (opt, default false)
%   keepUnbalanced      include reactions cannot be balanced, usually
%                       because they are polymeric reactions or with
%                       specific difficulty in balancing class structures
%                       (opt, default false)
%   keepUndetermined    include reactions that have substrates lack chemical
%                       structures or with non-numerical coefficients (e.g. n+1)
%                       (opt, default false)
%
%   Output:
%   model     a model structure generated from the database. The following
%             fields are filled
%             id:             'MetaCyc'
%             description:    'Automatically generated from MetaCyc database'
%             rxns:           MetaCyc reaction ids
%             rxnNames:       Name for each reaction entry
%             mets:           MetaCyc compound ids. If the equations use
%                             stoichiometry such as ID(n+1) then the whole
%                             expression is saved as the id
%             pwys:           MetaCyc pathway id
%             eccodes:        Corresponding ec-number if available
%             rxnMiriams:     Contains reaction specific information such as
%                             associated RHEA and KEGG reaction ids
%             S:              Stoichiometric matrix
%             lb:             -1000 for all reactions
%             ub:             1000 for all reactions
%             rev:            1 for reversible and 0 for irreversible. For
%                             reactions present in pathway maps the reversibility
%                             is taken from there
%             b:              0 for all metabolites
%             version:        MetaCyc database version
%
%   If the file metaCycRxns.mat is in the RAVEN/external/metacyc directory, it
%   will be directly loaded instead of parsing the MetaCyc data files and
%   pre-prepared lists of MetaCyc transport and undetermined reactions.
%
%   Usage: model=getRxnsFromMetaCyc(metacycPath,keepTransportRxns,keepUnbalanced,keepUndetermined)

%NOTE: This is how one entry looks in the file

% UNIQUE-ID - 4-HYDROXY-2-KETOPIMELATE-LYSIS-RXN
% TYPES - Small-Molecule-Reactions
% TYPES - Chemical-Reactions
% ATOM-MAPPINGS - (:NO-HYDROGEN-ENCODING (6 7 0 8 1 9 2 10 3 12 11 5 4) (((PYRUVATE 0 5) (SUCC-S-ALD 6 12)) ((CPD-804 0 12))))
% DBLINKS - (LIGAND-RXN "R01647" NIL |kr| 3650824949 NIL NIL)
% DBLINKS - (LIGAND-RXN "R01645" NIL |taltman| 3459474590 NIL NIL)
% EC-NUMBER - EC-4.1.2.52
% ^OFFICIAL? - T
% ENZYMATIC-REACTION - ENZRXN-3009
% IN-PATHWAY - 3-HYDROXYPHENYLACETATE-DEGRADATION-PWY
% LEFT - CPD-804
% ORPHAN? - :NO
% PHYSIOLOGICALLY-RELEVANT? - T
% REACTION-DIRECTION - LEFT-TO-RIGHT
% RIGHT - SUCC-S-ALD
% RIGHT - PYRUVATE
% //

% A line that contains only '//' separates each object.

if nargin<2
    keepTransportRxns=false;
end
if nargin<3
    keepUnbalanced=false;
end
if nargin<4
    keepUndetermined=false;
end

%Check if the reactions have been parsed before and saved. Directly load
%the model if so.
[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
rxnsFile=fullfile(ravenPath,'external','metacyc','metaCycRxns.mat');
metaCycRxnFile='reactions.dat';
metaCycPwyFile='pathway-links.dat';

if exist(rxnsFile, 'file')
    fprintf(['Importing MetaCyc reactions from ' strrep(rxnsFile,'\','/') '... ']);
    load(rxnsFile);
    fprintf('done\n');
else
    fprintf(['Cannot locate ' strrep(rxnsFile,'\','/') '\nNow try to generate it from local MetaCyc data files...\n']);
    if ~exist(fullfile(metacycPath,metaCycRxnFile),'file') || ~exist(fullfile(metacycPath,metaCycPwyFile),'file')
        EM=fprintf(['The files of reactions or pathways cannot be located, and should be downloaded from MetaCyc.\n']);
        dispEM(EM);
    else
        metaCycRxns.id='MetaCyc';
        metaCycRxns.description='Automatically generated from MetaCyc database';
        
        %Get pathway names and add them to the field of subSystems
        fid = fopen(fullfile(metacycPath,metaCycPwyFile), 'r');
        
        %Keeps track of how many pathways that have been added
        pwyCounter=0;
        pwys=cell(10000,1);
        pwyNames=cell(10000,1);
        
        %Loop through the file
        while 1
            %Get the next line
            tline = fgetl(fid);
            
            %Abort at end of file
            if ~ischar(tline)
                break;
            end
            
            %Read in pathway ids and names
            if ~strcmp(tline(1),'#')
                pwyCounter=pwyCounter+1;
                vars = regexp(tline, '\t', 'split');
                pwys{pwyCounter}=vars{1};
                pwyNames{pwyCounter}=vars{2};
                
                %Romve HTML symbols
                pwyNames{pwyCounter}=removeHTMLcodes(pwyNames{pwyCounter});
                
            end
        end
        fclose(fid);
        pwys=pwys(1:pwyCounter);
        pwyNames=pwyNames(1:pwyCounter);
        
        %Preallocate memory for 50000 reactions
        metaCycRxns.rxns=cell(50000,1);
        metaCycRxns.rxnNames=cell(50000,1);
        metaCycRxns.eccodes=cell(50000,1);
        metaCycRxns.subSystems=cell(50000,1);
        metaCycRxns.pwys=cell(50000,1);
        metaCycRxns.rxnMiriams=cell(50000,1);
        metaCycRxns.rxnReferences=cell(50000,1);
        metaCycRxns.rev=ones(50000,1); %reversibility;
        
        rxnLinks.metacyc=cell(10000,1);
        rxnLinks.kegg=cell(10000,1);
        rxnLinks.check=cell(10000,1);
        isSpontaneous=false(10000,1); %spontaneous;
        UNBALANCED=false(10000,1);
        UNDETERMINED=false(10000,1);
        TRANSPORT={};                 %transport reactions;
        
        metaCycRxns.equations=cell(50000,1); %reaction equations
        left=cell(50000,1); %Temporarily stores the equations
        right=cell(50000,1); %Temporarily stores the equations
        
        %First load information on reaction ID, reaction name,, pathway,
        %and ec-number
        fid = fopen(fullfile(metacycPath,metaCycRxnFile), 'r');
        
        %Keeps track of how many reactions that have been added
        rxnCounter=0;
        addSpont=false;
        dbLinkCounter=0;
        
        %Loop through the file
        while 1
            %Get the next line
            tline = fgetl(fid);
            
            %Abort at end of file
            if ~ischar(tline)
                break;
            end

            % Get the version of MetaCyc database
            if numel(tline)>11 && strcmp(tline(1:11),'# Version: ')
                version=tline(12:end);
            end

            %Check if it is a new reaction
            if numel(tline)>12 && strcmp(tline(1:12),'UNIQUE-ID - ')
                rxnCounter=rxnCounter+1;
                nPwys=0;
                
                %Add empty strings where there should be such
                metaCycRxns.rxnNames{rxnCounter}='';
                metaCycRxns.eccodes{rxnCounter}='';
                metaCycRxns.subSystems{rxnCounter}='';
                metaCycRxns.pwys{rxnCounter}='';
                metaCycRxns.equations{rxnCounter}='';
                metaCycRxns.rxnReferences{rxnCounter}='';
                reverse=0;   %Value for reversing the equation when necessary
                
                % A complex evaluation system for generating the equations
                left{rxnCounter}='';
                right{rxnCounter}='';
                coefficient='';
                templeft='';
                tempright='';
                % might be simplified with a better algorithem
                
                %Add reaction ID
                metaCycRxns.rxns{rxnCounter}=tline(13:end);
                
            end
            
            %Add name
            if numel(tline)>14 && strcmp(tline(1:14),'COMMON-NAME - ')
                metaCycRxns.rxnNames{rxnCounter}=tline(15:end);
                
                %Romve HTML symbols
                metaCycRxns.rxnNames{rxnCounter}=removeHTMLcodes(metaCycRxns.rxnNames{rxnCounter});
            end
            
            %Add eccodes
            if numel(tline)>15 && strcmp(tline(1:15),'EC-NUMBER - EC-')
                if isempty(metaCycRxns.eccodes{rxnCounter})
                    metaCycRxns.eccodes{rxnCounter}=strcat('ec-code/',tline(16:end));
                else
                    metaCycRxns.eccodes{rxnCounter}=strcat(metaCycRxns.eccodes{rxnCounter},';ec-code/',tline(16:end));
                end
            end
            if numel(tline)>16 && strcmp(tline(1:16),'EC-NUMBER - |EC-')
                if isempty(metaCycRxns.eccodes{rxnCounter})
                    metaCycRxns.eccodes{rxnCounter}=strcat('ec-code/',tline(17:end-1));
                else
                    metaCycRxns.eccodes{rxnCounter}=strcat(metaCycRxns.eccodes{rxnCounter},';ec-code/',tline(17:end-1));
                end
            end
            
            %Add pathway id and assign pathway name to subSystem field
            if numel(tline)>13 && strcmp(tline(1:13),'IN-PATHWAY - ')
                if isempty(metaCycRxns.pwys{rxnCounter})
                    metaCycRxns.pwys{rxnCounter}=tline(14:end);
                else
                    metaCycRxns.pwys{rxnCounter}=strcat(metaCycRxns.pwys{rxnCounter},';',tline(14:end));
                end
                
                [x, y]=ismember(tline(14:end),pwys);
                if x
                    metaCycRxns.subSystems{rxnCounter,1}{1,numel(metaCycRxns.subSystems{rxnCounter,1})+1}=pwyNames{y};
                end
            end
            
            %Add references (pubmed ids)
            if numel(tline)>12 && strcmp(tline(1:12),'CITATIONS - ')
                if isempty(metaCycRxns.rxnReferences{rxnCounter})
                    metaCycRxns.rxnReferences{rxnCounter}=strcat('pubmed/',tline(13:end));
                else
                    metaCycRxns.rxnReferences{rxnCounter}=strcat(metaCycRxns.rxnReferences{rxnCounter},';pubmed/',tline(13:end));
                end
            end
            
            %Add Miriam info (cross-links other database)
            if numel(tline)>11 && strcmp(tline(1:11),'DBLINKS - (')
                dblink=tline(12:end);
                %KEGG reaction id
                if strcmp(dblink(1:12),'LIGAND-RXN "')
                    dblink=dblink(13:end);
                    s=strfind(dblink,'"');
                    if any(s)
                        dblink=dblink(1:s-1);
                    end
                    
                    if isstruct(metaCycRxns.rxnMiriams{rxnCounter})
                        addToIndex=numel(metaCycRxns.rxnMiriams{rxnCounter}.name)+1;
                    else
                        addToIndex=1;
                    end
                    tempStruct=metaCycRxns.rxnMiriams{rxnCounter};
                    tempStruct.name{addToIndex,1}='kegg.reaction';
                    tempStruct.value{addToIndex,1}=dblink;
                    metaCycRxns.rxnMiriams{rxnCounter}=tempStruct;
                    
                    %For generating the rxnLinks structure
                    dbLinkCounter=dbLinkCounter+1;
                    rxnLinks.metacyc{dbLinkCounter}=metaCycRxns.rxns{rxnCounter};
                    rxnLinks.kegg{dbLinkCounter}=dblink;
                    rxnLinks.check{dbLinkCounter}=strcat(metaCycRxns.rxns{rxnCounter},dblink);
                end
                
                %RHEA reaction id
                if strcmp(dblink(1:6),'RHEA "')
                    dblink=dblink(7:end);
                    s=strfind(dblink,'"');
                    if any(s)
                        dblink=dblink(1:s-1);
                    end
                    
                    if isstruct(metaCycRxns.rxnMiriams{rxnCounter})
                        addToIndex=numel(metaCycRxns.rxnMiriams{rxnCounter}.name)+1;
                    else
                        addToIndex=1;
                    end
                    tempStruct=metaCycRxns.rxnMiriams{rxnCounter};
                    tempStruct.name{addToIndex,1}='rhea';
                    tempStruct.value{addToIndex,1}=dblink;
                    metaCycRxns.rxnMiriams{rxnCounter}=tempStruct;
                end
            end
            
            if numel(tline)>21 && strcmp(tline(1:21),'REACTION-DIRECTION - ')
                rxnDirection=tline(22:end);
                switch(rxnDirection)
                    case 'IRREVERSIBLE-LEFT-TO-RIGHT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                    case 'LEFT-TO-RIGHT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                    case 'PHYSIOL-LEFT-TO-RIGHT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                    case 'IRREVERSIBLE-RIGHT-TO-LEFT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                        reverse=1;
                    case 'RIGHT-TO-LEFT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                        reverse=1;
                    case 'PHYSIOL-RIGHT-TO-LEFT'
                        metaCycRxns.rev(rxnCounter,1)=0;
                        reverse=1;
                end
            end
            
            %Tag transport reactions
            if strcmp(tline,'TYPES - Transport-Reactions')
                TRANSPORT=[TRANSPORT;metaCycRxns.rxns{rxnCounter}];
            end
            
            %Add spontaneous
            if strcmp(tline,'SPONTANEOUS? - T')
                %metaCycRxns.spontaneous(rxnCounter,1)=1;
                isSpontaneous(rxnCounter)=true;
            end
            
            %Extract mass-balance status
            if numel(tline)>27 && strcmp(tline(1:27),'REACTION-BALANCE-STATUS - :')
                if isequal(tline(28:35), 'UNBALANC')
                    UNBALANCED(rxnCounter)=true;
                elseif isequal(tline(28:35), 'UNDETERM')
                    UNDETERMINED(rxnCounter)=true;
                end
            end
            
            %Add left side equation
            if numel(tline)>7 && strcmp(tline(1:7),'LEFT - ')
                if strcmp(left{rxnCounter},'')
                    if strcmp(templeft,'')
                        templeft=tline(8:end); % this is the real first, save to a temp variable
                    else % count in the coefficient here
                        if strcmp(coefficient,'')
                            left{rxnCounter}=templeft;
                        else
                            left{rxnCounter}=strcat(coefficient,32,templeft);
                            coefficient='';
                        end
                        templeft=tline(8:end);
                    end
                else
                    if strcmp(coefficient,'')
                        left{rxnCounter}=strcat(left{rxnCounter},' +',32,templeft);
                    else
                        left{rxnCounter}=strcat(left{rxnCounter},' +',32,coefficient,32,templeft);
                        coefficient='';
                    end
                    templeft=tline(8:end);
                end
            end
            
            %Add right side equation
            if numel(tline)>8 && strcmp(tline(1:8),'RIGHT - ')
                if strcmp(right{rxnCounter},'')
                    if strcmp(tempright,'')
                        
                        % a complicated process for the last left
                        % metabolite
                        if strcmp(coefficient,'')
                            if strcmp(left{rxnCounter},'')
                                left{rxnCounter}=templeft;
                            else
                                left{rxnCounter}=strcat(left{rxnCounter},' +',32,templeft);
                            end
                            
                        else
                            if strcmp(left{rxnCounter},'')
                                left{rxnCounter}=strcat(coefficient,32,templeft);
                            else
                                left{rxnCounter}=strcat(left{rxnCounter},' +',32,coefficient,32,templeft);
                            end
                            coefficient='';
                        end
                        % process end
                        
                        tempright=tline(9:end); % this is the real first, save to a temp variable
                    else
                        if strcmp(coefficient,'')
                            right{rxnCounter}=tempright;
                        else
                            right{rxnCounter}=strcat(coefficient,32,tempright);
                            coefficient='';
                        end
                        tempright=tline(9:end);
                    end
                else
                    if strcmp(coefficient,'')
                        right{rxnCounter}=strcat(right{rxnCounter},' +',32,tempright);
                    else
                        right{rxnCounter}=strcat(right{rxnCounter},' +',32,coefficient,32,tempright);
                        coefficient='';
                    end
                    tempright=tline(9:end);
                end
                
            end
            
            if numel(tline)>15 && strcmp(tline(1:15),'^COEFFICIENT - ')
                coefficient=tline(16:end);
            end
            
            %Generate equation at the end of each object section
            if strcmp(tline,'//')
                
                % a complicated process for the last right metabolite,
                % sub-function is needed here
                if strcmp(coefficient,'')
                    if strcmp(right{rxnCounter},'')
                        right{rxnCounter}=tempright;
                    else
                        right{rxnCounter}=strcat(right{rxnCounter},' +',32,tempright);
                    end
                    
                else
                    if strcmp(right{rxnCounter},'')
                        right{rxnCounter}=strcat(coefficient,32,tempright);
                    else
                        right{rxnCounter}=strcat(right{rxnCounter},' +',32,coefficient,32,tempright);
                    end
                    coefficient='';
                end
                % process end
                
                %get the right direction symbol
                if metaCycRxns.rev(rxnCounter,1)
                    symbol = ' <=>';
                else
                    symbol = ' =>';
                end
                
                if reverse
                    metaCycRxns.equations{rxnCounter}=strcat(right{rxnCounter},symbol,32,left{rxnCounter});
                else
                    metaCycRxns.equations{rxnCounter}=strcat(left{rxnCounter},symbol,32,right{rxnCounter});
                end
                
                %Final equation check
                if strcmp(left{rxnCounter},'') || strcmp(right{rxnCounter},'')
                    rxnCounter=rxnCounter-1;
                    
                end
                
            end
            
        end
        %Close the file
        fclose(fid);
        
        %===
        UNBALANCED=metaCycRxns.rxns(UNBALANCED);
        UNDETERMINED=metaCycRxns.rxns(UNDETERMINED);
        isSpontaneous=metaCycRxns.rxns(isSpontaneous);
        
        %If too much space was allocated, shrink the model
        metaCycRxns.rxns=metaCycRxns.rxns(1:rxnCounter);
        metaCycRxns.rxnNames=metaCycRxns.rxnNames(1:rxnCounter);
        metaCycRxns.eccodes=metaCycRxns.eccodes(1:rxnCounter);
        metaCycRxns.equations=metaCycRxns.equations(1:rxnCounter);
        metaCycRxns.rxnMiriams=metaCycRxns.rxnMiriams(1:rxnCounter);
        metaCycRxns.rxnReferences=metaCycRxns.rxnReferences(1:rxnCounter);
        metaCycRxns.subSystems=metaCycRxns.subSystems(1:rxnCounter);
        metaCycRxns.pwys=metaCycRxns.pwys(1:rxnCounter);
        metaCycRxns.rev=metaCycRxns.rev(1:rxnCounter,:);
        
        rxnLinks.kegg=rxnLinks.kegg(1:dbLinkCounter);
        rxnLinks.metacyc=rxnLinks.metacyc(1:dbLinkCounter);
        rxnLinks.check=rxnLinks.check(1:dbLinkCounter);
        [~,index]=unique(rxnLinks.check);
        rxnLinks.kegg=rxnLinks.kegg(index);
        rxnLinks.metacyc=rxnLinks.metacyc(index);
        rxnLinks=rmfield(rxnLinks,'check');
        
        %Construct the S matrix and list of metabolites
        [S, mets, badRxns]=constructS(metaCycRxns.equations);
        metaCycRxns.S=S;
        metaCycRxns.mets=mets;
        
        %Add some stuff to get a correct model structure
        metaCycRxns.ub=ones(rxnCounter,1)*1000;
        metaCycRxns.lb=metaCycRxns.rev*-1000;
        metaCycRxns.c=zeros(rxnCounter,1);
        metaCycRxns.b=zeros(numel(metaCycRxns.mets),1);
        metaCycRxns.version=version;
        
        %Save the model structure
        save(rxnsFile,'metaCycRxns','rxnLinks','TRANSPORT','UNBALANCED','UNDETERMINED','isSpontaneous');
        fprintf(['New metaCycRxns.mat has been successfully updated!\n\n']);
    end
end

%Deal with reactions that are labeled as "TRANSPORT", "UNBALANCED", or
%"UNDETERMINED" (depending on settings).
model=metaCycRxns;
if keepTransportRxns==false
    model=removeReactions(model,intersect(TRANSPORT,model.rxns),true,true);
end
if keepUnbalanced==false
    model=removeReactions(model,intersect(UNBALANCED,model.rxns),true,true);
end
if keepUndetermined==false
    model=removeReactions(model,intersect(UNDETERMINED,model.rxns),true,true);
end
end

%Sub function for romving HTML symbols from the names of reaction and
%pathway
function newString=removeHTMLcodes(string)
string=regexprep(string,'<(\w+)>','');
string=regexprep(string,'</(\w+)>','');
string=regexprep(string,'[&;]','');
newString=string;
end
