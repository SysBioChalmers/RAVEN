classdef tCuration < RavenTestCase
% tCuration  Tests for the table-driven curation functions in curation/.
%
%   curateModelFromTables reads its mets/genes/rxns input from tab-delimited
%   files; the tests below write those files to a temporary directory.

    methods (Test)

        function curateModelFromTablesKeepsSupersetReaction(testCase)
            % A model reaction that uses the candidate's metabolites *and*
            % others is a different reaction, not a duplicate: adding
            % A + B -> C to a model holding A + B -> C + H must leave the
            % existing reaction alone and add the candidate separately.
            model = tCuration.toyModel();
            [coeffsF, infoF] = tCuration.writeRxnTables(testCase, ...
                {'A','B','C'}, [-1 -1 1]);
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyNumElements(out.rxns, 2);
            % The existing reaction keeps its name and its fourth metabolite
            testCase.verifyEqual(out.rxnNames{1}, 'ATP hydrolysis with proton');
            testCase.verifyNumElements(find(out.S(:,1)), 4);
            % The candidate is added as a new reaction
            testCase.verifyEqual(out.rxnNames{2}, 'newRxn');
            testCase.verifyNumElements(find(out.S(:,2)), 3);
        end

        function curateModelFromTablesOverwritesExactDuplicate(testCase)
            % An exact stoichiometric match is still treated as the same
            % reaction and has its annotation overwritten.
            model = tCuration.toyModel();
            [coeffsF, infoF] = tCuration.writeRxnTables(testCase, ...
                {'A','B','C','H'}, [-1 -1 1 1]);
            evalc(['out = curateModelFromTables(model,''none'',''rxnsCoeffs'',' ...
                'coeffsF,''rxnsInfo'',infoF);']);

            testCase.verifyNumElements(out.rxns, 1);
            testCase.verifyEqual(out.rxnNames{1}, 'newRxn');
        end

    end

    methods (Static, Access = private)

        function model = toyModel()
            % Single reaction: A + B -> C + H, all in one compartment.
            model = struct();
            model.id        = 'toy';
            model.rxns      = {'R1'};
            model.rxnNames  = {'ATP hydrolysis with proton'};
            model.mets      = {'m1';'m2';'m3';'m4'};
            model.metNames  = {'A';'B';'C';'H'};
            model.comps     = {'c'};
            model.compNames = {'cytosol'};
            model.metComps  = [1;1;1;1];
            model.S         = sparse([-1;-1;1;1]);
            model.lb        = 0;
            model.ub        = 1000;
            model.rev       = 0;
            model.c         = 0;
            model.b         = zeros(4,1);
            model.genes     = {};
            model.grRules   = {''};
            model.rxnGeneMat= sparse(1,0);
            model.metFormulas = {'';'';'';''};
            model.subSystems  = {{''}};
            model.eccodes     = {''};
            model.rxnNotes    = {''};
            model.rxnReferences = {''};
            model.rxnConfidenceScores = 0;
        end

        function [coeffsF, infoF] = writeRxnTables(testCase, metNames, coeffs)
            % Write a one-reaction rxnsCoeffs/rxnsInfo pair for 'newRxn'.
            d = tempname; mkdir(d);
            testCase.addTeardown(@() rmdir(d,'s'));
            coeffsF = fullfile(d,'coeffs.tsv');
            infoF   = fullfile(d,'info.tsv');

            fid = fopen(coeffsF,'w');
            fprintf(fid,'rxnIdx\trxnNames\tmetNames\tcomps\tcoefficient\n');
            for i=1:numel(metNames)
                fprintf(fid,'1\tnewRxn\t%s\tc\t%d\n', metNames{i}, coeffs(i));
            end
            fclose(fid);

            fid = fopen(infoF,'w');
            fprintf(fid,['rxnIdx\trxnNames\tgrRules\tlb\tub\trev\tsubSystems\t' ...
                'eccodes\trxnNotes\trxnReferences\trxnConfidenceScores\n']);
            fprintf(fid,'1\tnewRxn\t\t0\t1000\t0\t\t\t\t\t\n');
            fclose(fid);
        end

    end
end
