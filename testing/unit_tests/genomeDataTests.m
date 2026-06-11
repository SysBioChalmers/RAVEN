%run this test case with the command
%results = runtests('genomeDataTests.m')
function tests = genomeDataTests
tests = functiontests(localfunctions);
end

function dataDir = fixtureDir()
dataDir = fullfile(fileparts(which(mfilename)), 'test_data', 'genomeData');
end

function getGeneData_eukaryoteTest(testCase)
% The eukaryote GFF3 has gene -> mRNA -> CDS; the protein accession must
% come from the CDS protein_id (NP_), not the mRNA transcript (NM_), and
% repeated CDS lines for one protein must be collapsed to a single row.
T = getGeneData(fullfile(fixtureDir(), 'euk.gff'));
verifyEqual(testCase, height(T), 2);
i1 = find(strcmp(T.locus_tag, 'LT1'), 1);
verifyEqual(testCase, T.old_locus_tag{i1}, 'OLD1');
verifyEqual(testCase, T.GeneID{i1}, '111');
verifyEqual(testCase, T.gene_name{i1}, 'AAA');
verifyEqual(testCase, T.GenBank_protein{i1}, 'NP_111.1');
i2 = find(strcmp(T.locus_tag, 'LT2'), 1);
verifyEqual(testCase, T.GenBank_protein{i2}, 'NP_222.1');
end

function getGeneData_prokaryoteTest(testCase)
% The prokaryote GFF3 has gene -> CDS directly.
T = getGeneData(fullfile(fixtureDir(), 'prok.gff'));
verifyEqual(testCase, height(T), 1);
verifyEqual(testCase, T.locus_tag{1}, 'b0001');
verifyEqual(testCase, T.GeneID{1}, '333');
verifyEqual(testCase, T.gene_name{1}, 'ccc');
verifyEqual(testCase, T.GenBank_protein{1}, 'WP_333.1');
end

function getGeneData_noFileWrittenByDefaultTest(testCase)
% Without an outputFile, no .tsv is written.
prevFiles = dir('geneData.tsv');
getGeneData(fullfile(fixtureDir(), 'prok.gff'));
verifyEqual(testCase, numel(dir('geneData.tsv')), numel(prevFiles));
end

function processProteinFastaFile_renamesHeadersTest(testCase)
T   = getGeneData(fullfile(fixtureDir(), 'euk.gff'));
out = processProteinFastaFile(fullfile(fixtureDir(), 'euk_protein.faa'), T, 'locus_tag', tempdir);
c = onCleanup(@() delete(out));
headers = {};
fid = fopen(out, 'r');
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(line) && line(1) == '>'
        headers{end+1} = strtrim(line(2:end)); %#ok<AGROW>
    end
end
fclose(fid);
% NP_111.1 -> LT1, NP_222.1 -> LT2, NP_999.9 has no match and is kept
verifyTrue(testCase, ismember('LT1', headers));
verifyTrue(testCase, ismember('LT2', headers));
verifyTrue(testCase, ismember('NP_999.9 unmatched', headers));
end

function renameModelGenes_renamesGenesAndGrRulesTest(testCase)
model.id          = 'test';
model.rxns        = {'r1'; 'r2'};
model.genes       = {'LT1'; 'LT2'};
model.grRules     = {'LT1'; 'LT1 or LT2'};
model.rxnGeneMat  = sparse([1 0; 1 1]);

T = getGeneData(fullfile(fixtureDir(), 'euk.gff'));
evalc('m2 = renameModelGenes(model, T, ''locus_tag'', ''gene_name'');');

verifyEqual(testCase, sort(m2.genes), {'AAA'; 'BBB'});
verifyFalse(testCase, any(contains(m2.grRules, 'LT1')));
verifyTrue(testCase, any(contains(m2.grRules, 'AAA')));
verifyTrue(testCase, any(contains(m2.grRules, 'BBB')));
% rxnGeneMat rebuilt to the new gene count
verifyEqual(testCase, size(m2.rxnGeneMat, 2), numel(m2.genes));
end

function downloadGenomeData_onlineTest(testCase)
% Network test: download a small genome from NCBI Datasets. Skipped when
% there is no internet connection.
d = tempname; mkdir(d);
c = onCleanup(@() rmdir(d, 's'));
try
    [gff, faa] = downloadGenomeData('GCF_000005845.2', d, false);
catch ME
    assumeFail(testCase, ['Skipped, likely no internet connection: ' ME.message]);
    return
end
verifyTrue(testCase, isfile(gff));
verifyTrue(testCase, isfile(faa));
end
