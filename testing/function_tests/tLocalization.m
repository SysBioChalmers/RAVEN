classdef tLocalization < RavenTestCase
% tLocalization  Tests for the subcellular-localization score parsers in localization/.
%
%   Exercises parseScores for each supported source (DeepLoc 2, MULocDeep, the COMPARTMENTS
%   database and UniProt) plus defaultCompartmentMap, using small fixture files written to a
%   temporary location. No external predictor or network access is needed.

    methods (Test)

        function deeploc(testCase)
            f = i_tmp(testCase, '.csv');
            writelines([
                "Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Mitochondrion"
                "G1,Cytoplasm,,0.8,0.1,0.05"
                "G2,Mitochondrion,SP,0.05,0.15,0.9"], f);
            GSS = parseScores(f, 'predictor', 'deeploc');
            testCase.verifyEqual(sort(string(GSS.compartments(:)')), ...
                ["Cytoplasm" "Mitochondrion" "Nucleus"]);
            testCase.verifyEqual(i_score(GSS,'G1','Cytoplasm'), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(i_score(GSS,'G2','Mitochondrion'), 1, 'AbsTol', 1e-9);
        end

        function deeplocWithCompartmentMap(testCase)
            f = i_tmp(testCase, '.csv');
            writelines([
                "Protein_ID,Localizations,Signals,Cytoplasm,Nucleus,Mitochondrion,Plastid"
                "G1,Cytoplasm,,0.8,0.1,0.05,0.4"], f);
            GSS = parseScores(f, 'predictor', 'deeploc', 'compartmentMap', defaultCompartmentMap());
            % labels mapped to model ids; Plastid (no fungal equivalent) dropped
            testCase.verifyEqual(sort(string(GSS.compartments(:)')), ["c" "m" "n"]);
            testCase.verifyEqual(i_score(GSS,'G1','c'), 1, 'AbsTol', 1e-9);
        end

        function mulocdeep(testCase)
            f = i_tmp(testCase, '.tsv');
            t = char(9);
            writelines([
                string(['protein' t 'Cytoplasm' t 'Mitochondrion' t 'Peroxisome'])
                string(['G1' t '0.2' t '0.7' t '0.1'])
                string(['G2' t '0.9' t '0.05' t '0.05'])], f);
            GSS = parseScores(f, 'predictor', 'mulocdeep', 'compartmentMap', defaultCompartmentMap());
            testCase.verifyEqual(sort(string(GSS.compartments(:)')), ["c" "m" "p"]);
            testCase.verifyEqual(i_score(GSS,'G1','m'), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(i_score(GSS,'G2','c'), 1, 'AbsTol', 1e-9);
        end

        function compartments(testCase)
            f = i_tmp(testCase, '.tsv');
            t = char(9);
            writelines([
                string(['YGL001C' t 'ERG26' t 'GO:0005739' t 'Mitochondrion' t '4.5'])
                string(['YGL001C' t 'ERG26' t 'GO:0005829' t 'Cytosol'       t '2.0'])
                string(['YGL002W' t 'ABC1'  t 'GO:0005777' t 'Peroxisome'    t '3.0'])
                string(['YGL002W' t 'ABC1'  t 'GO:0099999' t 'Plastid'       t '5.0'])], f);
            GSS = parseScores(f, 'predictor', 'compartments', ...
                              'compartmentMap', defaultCompartmentMap());
            testCase.verifyEqual(sort(string(GSS.compartments(:)')), ["c" "m" "p"]);
            testCase.verifyEqual(i_score(GSS,'YGL001C','m'), 1, 'AbsTol', 1e-9);          % 4.5 → 1
            testCase.verifyEqual(i_score(GSS,'YGL001C','c'), 2.0/4.5, 'AbsTol', 1e-9);
            % minConfidence drops the weak cytosol annotation
            GSS2 = parseScores(f, 'predictor', 'compartments', ...
                               'compartmentMap', defaultCompartmentMap(), 'minConfidence', 3.0);
            testCase.verifyFalse(any(strcmp(GSS2.compartments,'c')));
        end

        function uniprot(testCase)
            f = i_tmp(testCase, '.tsv');
            t = char(9);
            writelines([
                string(['Entry' t 'Gene Names (primary)' t 'Gene Names (ordered locus)' t ...
                        'Subcellular location [CC]'])
                string(['P00890' t 'CIT1' t 'YNR001C' t ...
                        'SUBCELLULAR LOCATION: Mitochondrion matrix {ECO:1}.'])
                string(['P12345' t 'GENEX' t 'YAL001C' t ...
                        'SUBCELLULAR LOCATION: Cytoplasm. Nucleus {ECO:2}. ' ...
                        'Note=Shuttles to the mitochondrion under stress.'])], f);
            GSS = parseScores(f, 'predictor', 'uniprot', 'idColumn', 'Gene Names (ordered locus)');
            testCase.verifyEqual(i_score(GSS,'YNR001C','m'), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(i_score(GSS,'YAL001C','c'), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(i_score(GSS,'YAL001C','n'), 1, 'AbsTol', 1e-9);
            % the "mitochondrion" mention inside Note=… is ignored
            testCase.verifyEqual(i_score(GSS,'YAL001C','m'), 0, 'AbsTol', 1e-9);
        end

        function uniprotWholeWordMatch(testCase)
            % 'Cytoplasmic vesicle' must NOT be mis-read as cytosol ('c') via substring match
            f = i_tmp(testCase, '.tsv');
            t = char(9);
            writelines([
                string(['Entry' t 'Subcellular location [CC]'])
                string(['YAA' t 'SUBCELLULAR LOCATION: Cytoplasmic vesicle.'])
                string(['YBB' t 'SUBCELLULAR LOCATION: Cytoplasm.'])], f);
            GSS = parseScores(f, 'predictor', 'uniprot');
            testCase.verifyEqual(i_score(GSS,'YBB','c'), 1, 'AbsTol', 1e-9);
            testCase.verifyFalse(any(strcmp(GSS.genes,'YAA')));   % no real location -> dropped
        end

        function uniprotMultiIsoformNote(testCase)
            % a Note= between isoform blocks must not swallow the later isoform's location
            f = i_tmp(testCase, '.tsv');
            t = char(9);
            writelines([
                string(['Entry' t 'Subcellular location [CC]'])
                string(['YCC' t 'SUBCELLULAR LOCATION: [Isoform 1]: Mitochondrion {ECO:1}. ' ...
                        'Note=Regulated under stress. [Isoform 2]: Cytoplasm {ECO:2}.'])], f);
            GSS = parseScores(f, 'predictor', 'uniprot');
            testCase.verifyEqual(i_score(GSS,'YCC','m'), 1, 'AbsTol', 1e-9);
            testCase.verifyEqual(i_score(GSS,'YCC','c'), 1, 'AbsTol', 1e-9);   % isoform 2 kept
        end

        function wolfRemoved(testCase)
            f = i_tmp(testCase, '.txt');
            writelines("anything", f);
            testCase.verifyError(@() parseScores(f, 'predictor', 'wolf'), 'RAVEN:badInput');
        end

    end
end

% ------------------------------------------------------------------------- helpers

function f = i_tmp(testCase, ext)
f = [tempname ext];
testCase.addTeardown(@() i_delete(f));
end

function i_delete(f)
if isfile(f); delete(f); end
end

function s = i_score(GSS, gene, comp)
gi = strcmp(GSS.genes, gene);
ci = strcmp(GSS.compartments, comp);
if ~any(gi) || ~any(ci)
    s = 0;
else
    s = GSS.scores(gi, ci);
end
end
