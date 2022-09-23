function model=combineMetaCycKEGGModels(metacycModel,keggModel)
% combineMetaCycKEGGModels
%	Combine MetaCyc and KEGG draft models into one model structure.
%
%	Input:
%	metacycModel	the reconstructed model from MetaCyc
%	keggModel       the reconstructed model from KEGG
%
%	Output:
%	model           a model structure generated by integrating information
%                   from draft models reconstructed using MetaCyc and KEGG
%                   databases. The 'rxnFrom/metFrom/geneFrom' fields are
%                   included to indicate the source.
%
%	Usage: model=combineMetaCycKEGGModels(metacycModel,keggModel)

%Just return the model
if nargin<2
    disp('Missing input model');
    return;
end

%Add MetaCyc model as template
model=metacycModel;
model.id='COMBINED';
model.name='Combined model from MetaCyc and KEGG draft models';

%Use MetaCyc model as template
model.rxnFrom=cell(numel(model.rxns),1);
model.rxnFrom(:)={'MetaCyc'};
model.metFrom=cell(numel(model.mets),1);
model.metFrom(:)={'MetaCyc'};
model.geneFrom=cell(numel(model.genes),1);
model.geneFrom(:)={'MetaCyc'};
model.grRulesKEGG=cell(numel(model.rxns),1);
model.grRulesKEGG(:)={''};

%Remove the S and rxnGeneMat matrix
if isfield(model,'S')
    model=rmfield(model,'S');
end
if isfield(model,'rxnGeneMat')
    model=rmfield(model,'rxnGeneMat');
end

%Arrange the genes related fields

%Identify the shared genes between KEGG and MetaCyc
sharedGenes=intersect(keggModel.genes,model.genes);
[~, b]=ismember(sharedGenes,model.genes);
model.geneFrom(b)={'Both'};

%Add unique genes from KEGG and update geneFrom
C=setdiff(keggModel.genes,model.genes);
model.genes=[model.genes;C];
geneFrom=cell(numel(C),1);
geneFrom(:)={'KEGG'};
model.geneFrom=[model.geneFrom;geneFrom];

%Prepare for matching grRules
if isfield(keggModel,'grRules')
    keggModel.grRules=strrep(keggModel.grRules,'(','');
    keggModel.grRules=strrep(keggModel.grRules,')','');
end

%Replace rxns in KEGG model with the corresponding ones in MetaCyc by
%updating values rxns field to MetaCyc version whenever possible

%Collect the ones found in MetaCyc as mappedRxns
mappedRxns=[];
rxnsToMove=[];
grRulesToMove={}; %Storing the grRules to be moved to MetaCyc model

%Read in dbLinks, note that linkMetaCycKEGGRxns need to be run in advance
ravenPath=findRAVENroot();
load(fullfile(ravenPath,'external','metacyc','metaCycRxns.mat'));
num=0;
numToMove=0;
%Loop through the KEGG model and replace the rxn id from KEGG to MetaCyc
%version
for i=1:numel(keggModel.rxns)
    [~, b]=ismember(rxnLinks.kegg,keggModel.rxns{i});
    I=find(b);
    
    if numel(I)==1
        %Find out the corresponding MetaCyc reactions
        num=num+1;
        %Record the mapped rxn index in KEGG model and remove them later
        mappedRxns=[mappedRxns,i];
        %keggModel.rxns{i}=rxnLinks.metacyc{I(1)};
        
        [c, d]=ismember(rxnLinks.metacyc{I(1)},model.rxns);
        if c
            model.rxnFrom{d}='Both';
            %Combine the grRule info Check if KEGG grRules equals or a
            %subset of MetaCyc grRules If not, save the KEGG grRules to a
            %new field grRulesKEGG for manual curation
            k=strfind(model.grRules{d},keggModel.grRules{i});
            if isempty(k)
                if isempty(model.grRulesKEGG{d})
                    model.grRulesKEGG{d}=keggModel.grRules{i};
                else
                    model.grRulesKEGG{d}=strcat(model.grRulesKEGG{d},{' or '},keggModel.grRules{i});
                end
            end
        else
            %For the KEGG rxns have MetaCyc version but not in MetaCyc
            %draft model
            [~, y]=ismember(rxnLinks.metacyc{I(1)},metaCycRxns.rxns);
            
            %Check if this reaction is repitetive and save the grRules
            [Repeat, Index]=ismember(y,rxnsToMove);
            if ~Repeat
                rxnsToMove=[rxnsToMove;y];  %Record rxns to be moved from KEGG to MetaCyc model
                if isempty(keggModel.grRules{i})
                    %Rxns may have empty grRules (e.g. spontaneous) that
                    %are identified and added back here
                    grRulesToMove=[grRulesToMove;{''}];
                else
                    grRulesToMove=[grRulesToMove;keggModel.grRules{i}];
                end
                numToMove=numToMove+1;
            else
                %If this reaction has been recorded, append the grRules if
                %it is different
                k=strfind(grRulesToMove{Index},keggModel.grRules{i});
                if isempty(k)
                    grRulesToMove{Index}=strcat(grRulesToMove{Index},{' or '},keggModel.grRules{i});
                end
            end
        end
    else
        %Here are the case of one KEGG rxn id link to several MetaCyc ids
        %for j=2:numel(I) Ignore this issue for now, and resolve it later,
        %because it never happened This case can be solved by adding above
        %section into a loop below disp(I(j)); end
    end
end
fprintf('NOTE: A total of %d reactions in the KEGG model were mapped to MetaCyc.\n',num);
fprintf('NOTE: %d reactions already in MetaCyc model, %d will be combined.\n',num-numToMove,numToMove);

%Append mapped MetaCyc rxns in KEGG model that are absent from MetaCyc
%model
rxnFrom=cell(numel(rxnsToMove),1);
rxnFrom(:)={'KEGG'};
model.rxnFrom=[model.rxnFrom;rxnFrom];
model.rxns=[model.rxns;metaCycRxns.rxns(rxnsToMove)];

rxnNames=metaCycRxns.rxnNames(rxnsToMove);
%Replace rxnNames with those from metaCycEnzymes
load('metaCycEnzymes.mat');
for i=1:numel(rxnsToMove)
    [a, b]=ismember(metaCycRxns.rxns{rxnsToMove(i)},metaCycEnzymes.rxns);
    if a
        rxnNames{i}=metaCycEnzymes.rxnNames{b};
    end
end
model.rxnNames=[model.rxnNames;rxnNames];

model.eccodes=[model.eccodes;metaCycRxns.eccodes(rxnsToMove)];
model.subSystems=[model.subSystems;metaCycRxns.subSystems(rxnsToMove)];
model.rxnMiriams=[model.rxnMiriams;metaCycRxns.rxnMiriams(rxnsToMove)];
model.rxnReferences=[model.rxnReferences;metaCycRxns.rxnReferences(rxnsToMove)];
model.rev=[model.rev;metaCycRxns.rev(rxnsToMove)];
model.equations=[model.equations;metaCycRxns.equations(rxnsToMove)];
model.lb=[model.lb;metaCycRxns.lb(rxnsToMove)];
model.ub=[model.ub;metaCycRxns.ub(rxnsToMove)];
model.c=[model.c;metaCycRxns.c(rxnsToMove)];
model.grRules=[model.grRules;grRulesToMove];
%model.grRules=[model.grRules;keggModel.grRules(rxnsToMove(:,2))];
%model.grRulesKEGG=[model.grRulesKEGG;keggModel.grRules(rxnsToMove(:,2))];
rxnFields=cell(numel(rxnsToMove),1);
rxnFields(:)={''};
%model.rxnNotes=[model.rxnNotes;rxnFields];
model.grRulesKEGG=[model.grRulesKEGG;rxnFields];

%Remove all mapped reactions from KEGG model
pureKeggModel=removeReactions(keggModel,mappedRxns,true,true);
fprintf('NOTE: A global KEGG model with %d reactions and %d metabolites was obtained.\n',numel(pureKeggModel.rxns),numel(pureKeggModel.mets));

%Replace mets information in KEGG model with the corresponding ones in
%MetaCyc This includes updating all metabolite-related fields, except S
%matrix.
load('metaCycMets.mat');

if ~isfield(pureKeggModel,'metCharges')
    pureKeggModel.metCharges=NaN(numel(pureKeggModel.mets),1);
end

num=0;
for i=1:numel(pureKeggModel.mets)
    [a, b]=ismember(pureKeggModel.mets{i},metaCycMets.keggid);
    if a
        %For now only replace the met ids
        num=num+1;
        pureKeggModel.mets{i}=metaCycMets.mets{b};
        %Record the replaced mets into a vector for later
    end
end
fprintf('NOTE: A total of %d metabolites from the global KEGG model will be re-mapped to MetaCyc.\n',num);

%Generate equations for pureKeggModel
equationStrings=constructEquations(pureKeggModel,'',false,false,false,true);

%Add the pure KEGG reactions to MetaCyc model
rxnFrom=cell(numel(pureKeggModel.rxns),1);
rxnFrom(:)={'KEGG'};
model.rxnFrom=[model.rxnFrom;rxnFrom];
model.rxns=[model.rxns;pureKeggModel.rxns];
model.equations=[model.equations;equationStrings];
model.rxnNames=[model.rxnNames;pureKeggModel.rxnNames];
model.eccodes=[model.eccodes;pureKeggModel.eccodes];
model.subSystems=[model.subSystems;pureKeggModel.subSystems];
model.rxnMiriams=[model.rxnMiriams;pureKeggModel.rxnMiriams];
model.rev=[model.rev;pureKeggModel.rev];
model.lb=[model.lb;pureKeggModel.lb];
model.ub=[model.ub;pureKeggModel.ub];
model.c=[model.c;pureKeggModel.c];
model.grRules=[model.grRules;pureKeggModel.grRules];
rxnFields=cell(numel(pureKeggModel.rxns),1);
rxnFields(:)={''};
model.grRulesKEGG=[model.grRulesKEGG;rxnFields];
model.rxnReferences=[model.rxnReferences;rxnFields];
%model.rxnNotes=[model.rxnNotes;rxnFields];

%It could also be that the metabolite and reaction names are empty for some
%reasons. In that case, use the ID instead
%I=cellfun(@isempty,model.metNames); model.metNames(I)=model.mets(I);
I=cellfun(@isempty,model.rxnNames);
model.rxnNames(I)=model.rxns(I);

%Generate S matrix and mets
[S, mets, ~]=constructS(model.equations);
model.S=S;
model.mets=mets;

%Remap metabolite information
if isfield(model,'metFrom')
    model=rmfield(model,'metFrom');
end
model.metFrom=cell(numel(model.mets),1);
model.metFrom(:)={''};

if isfield(model,'metNames')
    model=rmfield(model,'metNames');
end
model.metNames=cell(numel(model.mets),1);
model.metNames(:)={''};

if isfield(model,'metFormulas')
    model=rmfield(model,'metFormulas');
end
model.metFormulas=cell(numel(model.mets),1);
model.metFormulas(:)={''};

if ~isfield(model,'metCharges')
    model=rmfield(model,'metCharges');
end
model.metCharges=NaN(numel(model.mets),1);

if ~isfield(model,'inchis')
    model=rmfield(model,'inchis');
end
model.inchis=cell(numel(model.mets),1);
model.inchis(:)={''};

if ~isfield(model,'metMiriams')
    model=rmfield(model,'metMiriams');
end
model.metMiriams=cell(numel(model.mets),1);
model.metMiriams(:)={''};

keggMets=getMetsFromKEGG();
for i=1:numel(model.mets)
    [a, b]=ismember(model.mets{i},metaCycMets.mets);
    if a
        %Replace the met related info
        model.mets{i}=metaCycMets.mets{b};
        model.metNames{i}=metaCycMets.metNames{b};
        model.metFormulas{i}=metaCycMets.metFormulas{b};
        model.metMiriams{i}=metaCycMets.metMiriams{b};
        model.inchis{i}=metaCycMets.inchis{b};
        model.metCharges(i,1)=metaCycMets.metCharges(b,1);
        num=num+1;
    else
        [c, d]=ismember(model.mets{i},keggMets.mets);
        if c
            model.mets{i}=keggMets.mets{d};
            model.metNames{i}=keggMets.metNames{d};
            model.metFormulas{i}=keggMets.metFormulas{d};
            model.metMiriams{i}=keggMets.metMiriams{d};
            model.inchis{i}=keggMets.inchis{d};
            num=num+1;
        end
    end
    
end

%It could also be that the metabolite and reaction names are empty for some
%reasons. In that case, use the ID instead
I=cellfun(@isempty,model.metNames);
model.metNames(I)=model.mets(I);

%Add confidenceScores to model in a rough way
if isfield(model,'rxnConfidenceScores')
    model=rmfield(model,'rxnConfidenceScores');
end
model.rxnConfidenceScores=NaN(numel(model.rxns),1);
model.rxnConfidenceScores(:)=2;

%Put all metabolites in one compartment called 's' (for system). This is
%done just to be more compatible with the rest of the code
model.comps={'s'};
model.compNames={'System'};
model.metComps=ones(numel(model.mets),1);
model.b=zeros(numel(model.mets),1);

%Create the new rxnGene matrix, get back to this later
model.rxnGeneMat=sparse(numel(model.rxns),numel(model.genes));

end
