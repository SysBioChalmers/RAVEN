%run this test case with the command
%results = runtests('modelCurationTests.m')
function tests = modelCurationTests
tests = functiontests(localfunctions);
end

function addExchangeRxnsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=addExchangeRxns(model,'both','6pgl_c');

%Perform the curation manually for comparison
modelManual=model;
modelManual.rxns(end+1)={'EXC_BOTH_6pgl_c'};
modelManual.lb(end+1)=-1000;
modelManual.ub(end+1)=1000;
modelManual.c(end+1)=0;
modelManual.rev(end+1)=1;
modelManual.rxnNames(end+1)={'6-phospho-D-glucono-1,5-lactone exchange (BOTH)'};
modelManual.grRules(end+1)={''};
modelManual.eccodes(end+1)={''};
modelManual.rxnGeneMat(end+1,:)=sparse(zeros(1,numel(modelManual.genes)));
modelManual.S(:,end+1)=sparse(zeros(numel(modelManual.mets),1));
modelManual.S(5,end)=-1;

verifyEqual(testCase,modelNew,modelManual)
end

function addTransportTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

evalc('modelNew=addTransport(model,''c'',''e'',{''6-phospho-D-glucono-1,5-lactone''},false,false,''test_'');');

%Perform the curation manually for comparison
modelManual=model;
modelManual.rxns(end+1)={'test_0001'};
modelManual.lb(end+1)=0;
modelManual.ub(end+1)=Inf;
modelManual.c(end+1)=0;
modelManual.rev(end+1)=0;
modelManual.rxnNames(end+1)={'6-phospho-D-glucono-1,5-lactone transport, Cytoplasm-Extraorganism'};
modelManual.grRules(end+1)={''};
modelManual.eccodes(end+1)={''};
modelManual.rxnGeneMat(end+1,:)=sparse(zeros(1,numel(modelManual.genes)));

modelManual.S(:,end+1)=sparse(zeros(numel(modelManual.mets),1));
modelManual.S(end+1,:)=sparse(zeros(1,numel(modelManual.rxns)));
modelManual.S([5,73],end)=[-1,1];

modelManual.mets(end+1)={'m_0001'};
modelManual.metNames(end+1)={'6-phospho-D-glucono-1,5-lactone'};
modelManual.metMiriams(end+1)=modelManual.metMiriams(5);
modelManual.b(end+1)=[0];
modelManual.metFormulas(end+1)={'C6H9O9P'};
modelManual.metComps(end+1)=[2];

verifyEqual(testCase,modelNew,modelManual)
end

function addGenesRavenTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

genesToAdd.genes={'testgene1','testgene2'};
genesToAdd.geneShortNames={'short1','short2'};
genesToAdd.geneMiriams{1}=struct('name',{{'testDB1';'testDB2'}},'value',...
    {{'testValue1';'testValue2'}});
genesToAdd.geneMiriams{2}=struct('name',{{'testDB2';'testDB1'}},'value',...
    {{'testValue2';'testValue1'}});

modelNew=addGenesRaven(model,genesToAdd);

%Perform the curation manually for comparison
modelManual=model;
modelManual.geneMiriams(1:numel(modelManual.genes),1)={[]};
modelManual.genes(end+1:end+2)=genesToAdd.genes;
modelManual.geneShortNames(end+1:end+2)=genesToAdd.geneShortNames;
modelManual.rxnGeneMat(:,end+1:end+2)=sparse(zeros(numel(modelManual.rxns),2));
modelManual.geneMiriams(end+1:end+2)=genesToAdd.geneMiriams;
verifyEqual(testCase,modelNew,modelManual)
end

function removeGenesTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=removeGenes(model,'b1817',true,true,false);

%Perform the curation manually for comparison
modelManual=model;
modelManual.genes(54)=[];
modelManual.geneShortNames(54)=[];
modelManual.rxnGeneMat(:,54)=[];
modelManual.rxnGeneMat(45,:)=[];
modelManual.grRules(50)={'(b2415 and b1621 and b2417 and b2416) or (b2415 and b2417 and b2416 and b1101)'};
modelManual.grRules(45)=[];
modelManual.rxns(45)=[];
modelManual.rxnNames(45)=[];
modelManual.lb(45)=[];
modelManual.ub(45)=[];
modelManual.c(45)=[];
modelManual.rev(45)=[];
modelManual.S(:,45)=[];
modelManual.eccodes(45)=[];
verifyEqual(testCase,modelNew,modelManual)
end

function addMetsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

metsToAdd.metNames={'metaboliteName','3-Phospho-D-glycerate'};
metsToAdd.compartments={'c','e'};
metsToAdd.metCharges=[2,0];
metsToAdd.metNotes={'this is just a test','another one'};
metsToAdd.metMiriams{1}=struct('name',{{'testDB1';'testDB2'}},'value',...
    {{'testValue1';'testValue2'}});
metsToAdd.metMiriams{2}=struct('name',{{'testDB1'}},'value',...
    {{'testValue1'}});

evalc('modelNew=addMets(model,metsToAdd,true);');

%Perform the curation manually for comparison
modelManual=model;
modelManual.metNotes(1:numel(modelManual.mets),1)={''};
modelManual.metCharges(1:numel(modelManual.mets),1)=NaN(numel(modelManual.mets),1);
modelManual.mets(end+1:end+2)={'m_0001','m_0002'};
modelManual.metNames(end+1:end+2)=metsToAdd.metNames;
modelManual.metNotes(end+1:end+2)=metsToAdd.metNotes;
modelManual.metMiriams(end+1:end+2)=metsToAdd.metMiriams;
modelManual.metCharges(end+1:end+2)=[2,NaN];
modelManual.b(end+1:end+2)=[0,0];
modelManual.metFormulas(end+1:end+2)={'','C3H4O7P'};
modelManual.metComps(end+1:end+2)=[1,2];
modelManual.S(end+1:end+2,:)=sparse(zeros(2,numel(modelManual.rxns)));

verifyEqual(testCase,modelManual,modelNew)
end


function addMets_oneCompTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

metsToAdd.metNames={'metaboliteName','3-Phospho-D-glycerate'};
metsToAdd.compartments={'e'};
metsToAdd.metCharges=[2,0];
metsToAdd.metNotes={'this is just a test','another one'};
metsToAdd.metMiriams{1}=struct('name',{{'testDB1';'testDB2'}},'value',...
    {{'testValue1';'testValue2'}});
metsToAdd.metMiriams{2}=struct('name',{{'testDB1'}},'value',...
    {{'testValue1'}});

evalc('modelNew=addMets(model,metsToAdd,true);');

%Perform the curation manually for comparison
modelManual=model;
modelManual.metNotes(1:numel(modelManual.mets),1)={''};
modelManual.metCharges(1:numel(modelManual.mets),1)=NaN(numel(modelManual.mets),1);
modelManual.mets(end+1:end+2)={'m_0001','m_0002'};
modelManual.metNames(end+1:end+2)=metsToAdd.metNames;
modelManual.metNotes(end+1:end+2)=metsToAdd.metNotes;
modelManual.metMiriams(end+1:end+2)=metsToAdd.metMiriams;
modelManual.metCharges(end+1:end+2)=[2,NaN];
modelManual.b(end+1:end+2)=[0,0];
modelManual.metFormulas(end+1:end+2)={'','C3H4O7P'};
modelManual.metComps(end+1:end+2)=2;
modelManual.S(end+1:end+2,:)=sparse(zeros(2,numel(modelManual.rxns)));

verifyEqual(testCase,modelManual,modelNew)
end

function removeMetsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=removeMets(model,'Acetate',true,true,true,true);

%Perform the curation manually for comparison
modelManual=model;
modelManual.mets(6:7)=[];
modelManual.metNames(6:7)=[];
modelManual.metMiriams(6:7)=[];
modelManual.b(6:7)=[];
modelManual.metFormulas(6:7)=[];
modelManual.metComps(6:7)=[];
modelManual.S(6:7,:)=[];

modelManual.rxnGeneMat(20,:)=[];
modelManual.grRules(20)=[];
modelManual.rxns(20)=[];
modelManual.rxnNames(20)=[];
modelManual.lb(20)=[];
modelManual.ub(20)=[];
modelManual.c(20)=[];
modelManual.rev(20)=[];
modelManual.S(:,20)=[];
modelManual.eccodes(20)=[];
verifyEqual(testCase,modelNew,modelManual)
end

function addRxnsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');
rxnsToAdd.rxns='test1';
rxnsToAdd.equations='2-Oxoglutarate => TEST';
rxnsToAdd.rxnNames={'testRxn1'}; % test cell array
rxnsToAdd.rxnMiriams{1}=struct('name',{{'testDB1';'testDB2'}},'value',...
    {{'testValue1';'testValue2'}});
rxnsToAdd.grRules="b0008"; % test string
evalc('modelNew=addRxns(model,rxnsToAdd,2,''c'',true);');

rxnsToAdd.rxns='test2';
rxnsToAdd=rmfield(rxnsToAdd,'equations');
rxnsToAdd.mets={{'6pgl_c','atp_c'}};
rxnsToAdd.stoichCoeffs=[[-1,2]];
rxnsToAdd.lb=-1000;
rxnsToAdd.grRules='test';
evalc('modelNew=addRxns(modelNew,rxnsToAdd,1,''c'',true,true);');

%Perform the curation manually for comparison
modelManual=model;

modelManual.rxnMiriams=cell(numel(modelManual.rxns),1);
modelManual.rxns(end+1:end+2)={'test1';'test2'};
modelManual.mets(end+1)={'m_0001'};
modelManual.S(end+1,:)=zeros(1,numel(modelManual.rxns)-2);
modelManual.S(:,end+1:end+2)=zeros(numel(modelManual.mets),2);
modelManual.S(end,96)=1;
modelManual.S([14,73],end-1)=[-1,1];
modelManual.S([5,17],end)=[-1,2];
modelManual.lb(end+1:end+2)=[0;-1000];
modelManual.ub(end+1:end+2)=[Inf;Inf];
modelManual.rev(end+1:end+2)=[0;1];
modelManual.c(end+1:end+2)=[0;0];
modelManual.b(end+1)=[0];
modelManual.rxnNames(end+1:end+2)={'testRxn1';'testRxn1'};
modelManual.grRules(end+1:end+2)={'b0008';'test'};

modelManual.rxnGeneMat(end+1:end+2,:)=zeros(2,numel(modelManual.genes));
modelManual.rxnGeneMat(:,end+1)=zeros(numel(modelManual.rxns),1);
modelManual.rxnGeneMat(97,end)=1;
modelManual.rxnGeneMat(end-1,1)=1;
modelManual.rxnGeneMat(end,138)=1;

modelManual.eccodes(end+1:end+2)={'';''};
modelManual.genes(end+1)={'test'};
modelManual.geneShortNames(end+1)={''};
modelManual.metNames(end+1)={'TEST'};
modelManual.metComps(end+1)=[1];
modelManual.metFormulas(end+1)={''};
modelManual.metMiriams(end+1)={[]};
modelManual.rxnMiriams(end+1:end+2)={struct('name',{{'testDB1';'testDB2'}},'value',...
    {{'testValue1';'testValue2'}})};

verifyEqual(testCase,modelManual,modelNew)
end

function removeReactionsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=removeMets(model,'Acetate',true,true,true,true);

%Perform the curation manually for comparison
modelManual=model;
modelManual.mets(6:7)=[];
modelManual.metNames(6:7)=[];
modelManual.metMiriams(6:7)=[];
modelManual.b(6:7)=[];
modelManual.metFormulas(6:7)=[];
modelManual.metComps(6:7)=[];
modelManual.S(6:7,:)=[];

modelManual.rxnGeneMat(20,:)=[];
modelManual.grRules(20)=[];
modelManual.rxns(20)=[];
modelManual.rxnNames(20)=[];
modelManual.lb(20)=[];
modelManual.ub(20)=[];
modelManual.c(20)=[];
modelManual.rev(20)=[];
modelManual.S(:,20)=[];
modelManual.eccodes(20)=[];
verifyEqual(testCase,modelNew,modelManual)
end

function getMetsInCompTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

evalc('[testOut.I, testOut.metNames]=getMetsInComp(model,''e'');');

testCheck.I=[false;false;false;false;false;false;true;false;true;false;false;false;false;false;true;false;false;false;false;true;false;false;false;false;true;false;false;false;true;true;false;true;false;false;true;false;true;false;true;false;false;true;false;true;false;false;true;false;true;false;false;false;false;false;true;false;true;false;false;false;true;false;true;false;false;false;false;false;false;true;false;false];
testCheck.metNames={'Acetate';'Acetaldehyde';'2-Oxoglutarate';'CO2';'Ethanol';'Formate';'D-Fructose';'Fumarate';'D-Glucose';'L-Glutamine';'L-Glutamate';'H2O';'H+';'D-Lactate';'L-Malate';'Ammonium';'O2';'Phosphate';'Pyruvate';'Succinate'};
verifyEqual(testCase,testOut,testCheck)
end

function getRxnsInCompTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

evalc('[testOut.I, testOut.rxnNames]=getRxnsInComp(model,''e'');');

testCheck.I=[20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39];
testCheck.rxnNames={'Acetate exchange';'Acetaldehyde exchange';'2-Oxoglutarate exchange';'CO2 exchange';'Ethanol exchange';'Formate exchange';'D-Fructose exchange';'Fumarate exchange';'D-Glucose exchange';'L-Glutamine exchange';'L-Glutamate exchange';'H+ exchange';'H2O exchange';'D-lactate exchange';'L-Malate exchange';'Ammonia exchange';'O2 exchange';'Phosphate exchange';'Pyruvate exchange';'Succinate exchange'};
verifyEqual(testCase,testOut,testCheck)
end

function getAllRxnsFromGenesTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');
evalc('testOut=getAllRxnsFromGenes(model,''ACONTa'');');

testCheck={'ACONTa';'ACONTb'};
verifyEqual(testCase,testOut,testCheck)
end

function changeRxnsTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

evalc('modelNew=changeRxns(model,''ACKr'',''2-Oxoglutarate => TEST'',2,''c'',true);');

%Perform the curation manually for comparison
modelManual=model;

modelManual.mets(end+1)={'m_0001'};
modelManual.S(end+1,:)=zeros(1,numel(modelManual.rxns));
modelManual.S(:,3)=0;
modelManual.S([14,73],3)=[-1,1];

modelManual.metNames(end+1)={'TEST'};
modelManual.metComps(end+1)=[1];
modelManual.metFormulas(end+1)={''};
modelManual.metMiriams(end+1)={[]};
modelManual.b(end+1)=[0];
modelManual.rev(3)=[0];

verifyEqual(testCase,modelNew,modelManual)
end

function removeBadRxnsTest(testCase)
%Load the model in RAVEN format
sourceDir=fileparts(fileparts(fileparts(which(mfilename))));
load(fullfile(sourceDir,'testing','unit_tests','test_data','ecoli_textbook.mat'), 'model');

[~,exchRxns]=getExchangeRxns(model);
modelNew=model;

modelNew.lb(exchRxns)=[0];
modelNew.ub(exchRxns)=[0];
modelNew.lb(modelNew.lb>0)=0;
modelNew.S(17,:)=1;
[~, testOut]=removeBadRxns(modelNew);
evalc('[~, testOut]=removeBadRxns(modelNew,3,modelNew.metComps==2);');

import matlab.unittest.constraints.AnyCellOf;
import matlab.unittest.constraints.ContainsSubstring;
verifyThat(testCase,AnyCellOf({'SUCDi','FRD7'}),ContainsSubstring(testOut{:}),'test')
%Either SUCDi or FRD7 should be given as output
end

function changeGrRulesTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=changeGrRules(model,{'ACKr','ACONTa'},{'b2296 and b1849';'b1849'},true);
modelNew=changeGrRules(modelNew,'ACONTa','b1849',false);

%Perform the curation manually for comparison
modelManual=model;
modelManual.grRules(3:4)={'b2296 and b1849';'b1849 or b1849'};
modelManual.rxnGeneMat(3:4,:)=0;
modelManual.rxnGeneMat(3,[57,76])=1;
modelManual.rxnGeneMat(4,57)=1;

verifyEqual(testCase,modelNew,modelManual)
end

function copyToCompsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=copyToComps(model,{'p','e'},'ACKr');

%Perform the curation manually for comparison
modelManual=model;

modelManual.rxns(end+1:end+2)={'ACKr_p';'ACKr_e'};
modelManual.mets(end+1:end+7)={'ac_c_p';'actp_c_p';'adp_c_p';'atp_c_p';'actp_c_e';'adp_c_e';'atp_c_e'};
modelManual.S(end+1:end+7,:)=zeros(7,numel(modelManual.rxns)-2);
modelManual.S(:,end+1:end+2)=zeros(numel(modelManual.mets),2);

modelManual.S([73,74,75,76],end-1)=[-1,1,1,-1];
modelManual.S([7,77,78,79],end)=[-1,1,1,-1];

modelManual.lb(end+1:end+2)=[-1000;-1000];
modelManual.ub(end+1:end+2)=[1000;1000];
modelManual.rev(end+1:end+2)=[1;1];
modelManual.c(end+1:end+2)=[0;0];
modelManual.rxnNames(end+1:end+2)={'acetate kinase';'acetate kinase'};
modelManual.grRules(end+1:end+2)={'b2296 or b1849 or b3115';'b2296 or b1849 or b3115'};
modelManual.eccodes(end+1:end+2)={'2.7.2.1';'2.7.2.1'};

modelManual.rxnGeneMat(end+1:end+2,:)=zeros(2,numel(modelManual.genes));
modelManual.rxnGeneMat(end-1:end,[57,76,97])=1;

modelManual.metNames(end+1:end+7)={'Acetate';'Acetyl phosphate';'ADP';'ATP';'Acetyl phosphate';'ADP';'ATP'};
modelManual.metComps(end+1:end+7)=[3;3;3;3;2;2;2];
modelManual.metFormulas(end+1:end+7)={'C2H3O2';'C2H3O5P';'C10H12N5O10P2';'C10H12N5O13P3';'C2H3O5P';'C10H12N5O10P2';'C10H12N5O13P3'};
modelManual.metMiriams(end+1:end+7)=modelManual.metMiriams([6,12,13,17,12,13,17]);
modelManual.b(end+1:end+7)=[0];

modelManual.comps(end+1)={'p'};
modelManual.compNames(end+1)={'p'};

verifyEqual(testCase,modelNew,modelManual)
end

function setExchangeBoundsTest(testCase)
sourceDir = fileparts(which(mfilename));
load([sourceDir,'/test_data/ecoli_textbook.mat'], 'model');

modelNew=setExchangeBounds(model,{'ac_e';'akg_e'},-500,500);

%Perform the curation manually for comparison
modelManual=model;

exchRxns=[21;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39];
modelManual.lb(exchRxns)=0;
modelManual.lb([20,22])=-500;
modelManual.ub([20,22])=500;

verifyEqual(testCase,modelNew,modelManual)
end

function addRxnsGenesMetsTest(testCase)
sourceDir = fileparts(which(mfilename));
load(fullfile(sourceDir,'test_data','ecoli_textbook.mat'), 'model');

sbmlFile=fullfile(sourceDir,'..','..','tutorial','empty.xml');
evalc('modelEmpty=importModel(sbmlFile)'); % Repress warnings

evalc('modelNew=addRxnsGenesMets(model,modelEmpty,''r1'',true);');

%Perform the curation manually for comparison
modelManual=model;

modelManual.rxns(end+1)={'r1'};
modelManual.mets(end+1:end+3)={'m1';'m2';'m3'};
modelManual.S(end+1:end+3,:)=zeros(3,numel(modelManual.rxns)-1);
modelManual.S(:,end+1)=zeros(numel(modelManual.mets),1);
modelManual.S([42,73,74,75],end)=[-1,-1,1,1];

modelManual.lb(end+1)=[0];
modelManual.ub(end+1)=[1000];
modelManual.rev(end+1)=[0];
modelManual.c(end+1)=[0];
modelManual.rxnNames(end+1)={'Breakdown of sucrose (invertase)'};
modelManual.grRules(end+1)={'g1'};
modelManual.eccodes(end+1)={''};
modelManual.rxnNotes(1:numel(modelManual.rxns),1)={''};
modelManual.rxnNotes(end)={'Added via addRxnsGenesMets()'};
modelManual.rxnConfidenceScores(1:numel(modelManual.rxns),1)=NaN;
modelManual.rxnConfidenceScores(end)=[0];

modelManual.genes(end+1)={'g1'};
modelManual.geneShortNames(end+1)={''};

modelManual.rxnGeneMat(end+1,:)=zeros(1,numel(modelManual.genes)-1);
modelManual.rxnGeneMat(:,end+1)=zeros(numel(modelManual.rxns),1);
modelManual.rxnGeneMat(end,end)=1;

modelManual.metNames(end+1:end+3)={'sucrose';'glucose';'fructose'};
modelManual.metComps(end+1:end+3)=[2;2;2];
modelManual.metFormulas(end+1:end+3)={'C12H22O11';'C6H12O6';'C6H12O6'};
modelManual.metMiriams(end+1:end+3)={[];[];[]};
modelManual.b(end+1:end+3)=[0];

verifyEqual(testCase,modelNew,modelManual)
end
