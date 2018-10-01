function model=getMetsFromKEGG(keggPath)
% getMetsFromKEGG
%   Retrieves information on all metabolites stored in KEGG database
%
%   keggPath	if keggMets.mat is not in the RAVEN\external\kegg
%               directory, this function will attempt to read data from a
%               local FTP dump of the KEGG database. keggPath is the path
%               to the root of this database
%
%   model       a model structure generated from the database. The
%               following fields are filled
%   	id              'KEGG'
%   	description     'Automatically generated from KEGG database'
%   	mets            KEGG compound ids
%   	metNames        Compound name. Only the first name will be saved if
%                       there are several synonyms
%   	metMiriams      If there is a CHEBI id available, then that will be
%                       saved here
%   	inchis          InChI string for the metabolite
%   	metFormulas     The chemical composition of the metabolite. This
%                       will only be loaded if there is no InChI string
%
%   If the file keggMets.mat is in the RAVEN\external\kegg directory it
%   will be loaded instead of parsing of the KEGG files. If it does not
%   exist it will be saved after parsing of the KEGG files. In general, you
%   should remove the keggMets.mat file if you want to rebuild the model
%   structure from a newer version of KEGG.
%               
%   Usage: model=getMetsFromKEGG(keggPath)
%
%   Simonas Marcisauskas, 2018-07-25
%
%
% NOTE: This is how one entry looks in the file
%
% ENTRY       C00001                      Compound
% NAME        H2O;
%             Water
% FORMULA     H2O
% EXACT_MASS  18.0106
% MOL_WEIGHT  18.0153
% REMARK      Same as: D00001
% REACTION    R00001 R00002 R00004 R00005 R00009 R00010 R00011 R00017
%             R00022 R00024 R00025 R00026 R00028 R00036 R00041 R00044
%             (list truncated)
% ENZYME      1.1.1.1         1.1.1.22        1.1.1.23        1.1.1.115
%             1.1.1.132       1.1.1.136       1.1.1.170       1.1.1.186
%             (list truncated)
% BRITE       Therapeutic category of drugs in Japan [BR:br08301]
%             (list truncated)
% DBLINKS     CAS: 7732-18-5
%             PubChem: 3303
%             ChEBI: 15377
%             (list truncated)
%
% Then a lot of info about the positions of the atoms and so on. It is not
% certain that each metabolite follows this structure exactly.
%
% The file is not tab-delimited. Instead each label is 12 characters
% (except for '///').
%
% Check if the reactions have been parsed before and saved. If so, load the
% model.
%

if nargin<1
    keggPath='RAVEN/external/kegg';
end

[ST, I]=dbstack('-completenames');
ravenPath=fileparts(fileparts(fileparts(ST(I).file)));
metsFile=fullfile(ravenPath,'external','kegg','keggMets.mat');
if exist(metsFile, 'file')
    fprintf(['NOTE: Importing KEGG metabolites from ' strrep(metsFile,'\','/') '.\n']);
    load(metsFile);
else
    fprintf(['Cannot locate ' strrep(metsFile,'\','/') ' and will try to generate it from the local KEGG database.\n']);
    if ~exist(fullfile(keggPath,'compound'),'file') || ~exist(fullfile(keggPath,'compound.inchi'),'file')
        EM=fprintf(['The files ''compound'' and ''compound.inchi'' cannot be located at ' strrep(keggPath,'\','/') '/ and should be downloaded from the KEGG FTP.\n']);
        dispEM(EM);
    else
        %Add new functionality in the order specified in models
        model.id='KEGG';
        model.description='Automatically generated from KEGG database';
        
        %Preallocate memory for 30000 metabolites
        model.mets=cell(30000,1);
        model.metNames=cell(30000,1);
        model.metFormulas=cell(30000,1);
        model.metMiriams=cell(30000,1);
        
        %First load information on metabolite ID, metabolite name,
        %composition, and ChEBI
        
        fid = fopen(fullfile(keggPath,'compound'), 'r');
        
        %Keeps track of how many metabolites that have been added
        metCounter=0;
        
        %Loop through the file
        while 1
            %Get the next line
            tline = fgetl(fid);
            
            %Abort at end of file
            if ~ischar(tline)
                break;
            end
            
            %Skip '///'
            if numel(tline)<12
                continue;
            end
            
            %Check if it's a new reaction
            if strcmp(tline(1:12),'ENTRY       ')
                metCounter=metCounter+1;
                
                %Add empty strings where there should be such
                model.metNames{metCounter}='';
                model.metFormulas{metCounter}='';
                
                %Add compound ID (always 6 characters)
                model.mets{metCounter}=tline(13:18);
                
                %Add the KEGG id as metMiriams
                if length(model.mets{metCounter})==6
                    miriamStruct=model.metMiriams{metCounter};
                    if strcmp('G',model.mets{metCounter}(1))
                        miriamStruct.name{1,1}='kegg.glycan';
                    else
                        miriamStruct.name{1,1}='kegg.compound';
                    end
                    miriamStruct.value{1,1}=tline(13:18);
                    model.metMiriams{metCounter}=miriamStruct;
                end
            end
            
            %Add name
            if strcmp(tline(1:12),'NAME        ')
                %If there are synonyms, then the last character is ';'
                if strcmp(tline(end),';')
                    model.metNames{metCounter}=tline(13:end-1);
                    %Semicolon can also occur in the middle, separating
                    %several synonims in the same line
                    model.metNames{metCounter} = regexprep(model.metNames{metCounter},';.+','');
                elseif regexp(tline,';')
                    model.metNames{metCounter}=tline(13:end);
                    model.metNames{metCounter} = regexprep(model.metNames{metCounter},';.+','');
                else
                    model.metNames{metCounter}=tline(13:end);
                end
            end
            
            %Add composition
            if strcmp(tline(1:12),'FORMULA     ')
                model.metFormulas{metCounter}=tline(13:end);
            end
            
            %Add PubChem id
            if numel(tline)>21
                if strcmp(tline(13:21),'PubChem: ')
                    if isstruct(model.metMiriams{metCounter})
                        addToIndex=numel(model.metMiriams{metCounter}.name)+1;
                    else
                        addToIndex=1;
                    end
                    miriamStruct=model.metMiriams{metCounter};
                    miriamStruct.name{addToIndex,1}='pubchem.substance';
                    miriamStruct.value{addToIndex,1}=tline(22:end);
                    model.metMiriams{metCounter}=miriamStruct;
                end
            end
            
            %Add CHEBI id
            if numel(tline)>19
                if strcmp(tline(13:19),'ChEBI: ')
                    if isstruct(model.metMiriams{metCounter})
                        addToIndex=numel(model.metMiriams{metCounter}.name)+1;
                    else
                        addToIndex=1;
                    end
                    chebiIDs=strsplit(tline(20:end),' ');
                    miriamStruct=model.metMiriams{metCounter};
                    for i=1:numel(chebiIDs)
                        miriamStruct.name{addToIndex,1}='chebi';
                        miriamStruct.value{addToIndex,1}=strcat('CHEBI:',chebiIDs{i});
                        addToIndex=addToIndex+1;
                    end
                    model.metMiriams{metCounter}=miriamStruct;
                end
            end
        end
        
        %Close the file
        fclose(fid);
        
        %If too much space was allocated, shrink the model
        model.mets=model.mets(1:metCounter);
        model.metNames=model.metNames(1:metCounter);
        model.metFormulas=model.metFormulas(1:metCounter);
        model.metMiriams=model.metMiriams(1:metCounter);
        
        %Then load the InChI strings from another file. Not all metabolites
        %will be present in the list
        
        inchIDs=cell(numel(model.mets),1);
        inchis=cell(numel(model.mets),1);
        
        %The format is metID*tab*string
        
        fid = fopen(fullfile(keggPath,'compound.inchi'), 'r');
        
        %Loop through the file
        counter=1;
        while 1
            %Get the next line
            tline = fgetl(fid);
            
            %Abort at end of file
            if ~ischar(tline)
                break;
            end
            
            %Get the ID and the InChI
            inchIDs{counter}=tline(1:6);
            inchis{counter}=tline(14:end);
            counter=counter+1;
        end
        
        %Close the file
        fclose(fid);
        
        inchIDs=inchIDs(1:counter-1);
        inchis=inchis(1:counter-1);
        
        %Find the metabolites that had InChI strings and add them to the
        %model
        [a, b]=ismember(inchIDs,model.mets);
        
        %If there were mets with InChIs but that were not in the list
        if ~all(a)
            EM='Not all metabolites with InChI strings were found in the original list';
            disp(EM);
        end
        
        model.inchis=cell(numel(model.mets),1);
        model.inchis(:)={''};
        model.inchis(b)=inchis;
        
        %Remove composition if InChI was found
        model.metFormulas(b)={''};
        
        %Ensuring that all model.metMiriams.value consist only of strings,
        %no double
        for i=1:(numel(model.mets))
            for j=1:(numel(model.metMiriams{i}))
                if isa(model.metMiriams{i}.value{j},'double')
                    model.metMiriams{i}.value{j}=num2str(model.metMiriams{i}.value{j});
                end
            end
        end
        
        %Removing fronting and trailing whitespace from metNames
        model.metNames = deblank(model.metNames);
        
        %Fixing redundant metNames. The first occurence of particular
        %metabolite name is not changed, but starting from the second
        %occurence, original metabolite name is concatenated with KEGG
        %COMPOUND id between the brackets
        for i=1:(numel(model.metNames))
            if ~isempty(model.metNames{i})
                if sum(ismember(model.metNames(1:i-1),model.metNames(i)))>=1
                    model.metNames(i) = strcat(model.metNames(i), ' (', model.mets(i),')');
                end
            end
        end
        %Saves the model
        save(metsFile,'model');
    end
end
end
