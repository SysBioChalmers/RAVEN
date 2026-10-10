% tutorial7
%   This exercise shows how to use the biomass/ functions to calibrate a
%   model's biomass pseudoreaction against measured condition-specific
%   composition data: getBiomassFractions to see what the model currently
%   predicts, scaleBiomassFraction to rescale individual components to
%   measured values (with one component left to absorb whatever is left,
%   so the total stays at 1 g/gDW), and setGAM to apply a
%   growth-associated maintenance value. See Tutorial 7 on the RAVEN wiki
%   for more details: https://github.com/SysBioChalmers/RAVEN/wiki/Tutorials
%
%   A small hand-built model is used instead of a real GEM, so the
%   pseudoreaction layout is visible in full below. Real GEMs lay this out
%   the same way -- one pseudoreaction per macromolecular component,
%   feeding into a top-level biomass pseudoreaction -- and the file
%   biomassConditions_template.tsv in this folder is a template for the
%   condition data itself: one row per component, one column per named
%   condition, values in g/gDW. Real, published examples of the same
%   shape of data (per-metabolite abundances a modeller derives such
%   fractions from) are yeast-GEM's data/physiology/biomassComposition_*.tsv
%   files; the numbers used here are illustrative, not measured.

model = buildBiomassTemplateModel();
biomassConfig = buildBiomassTemplateConfig();

%What the template model's biomass pseudoreaction predicts before any
%calibration
fractions = getBiomassFractions(model, biomassConfig);
fprintf('Before calibration:\n');
disp(fractions);

%Load the condition data. Each column after "component" is a named
%condition; "carbohydrate" is deliberately not a row in the file -- it is
%the component left to balance the total back to 1 g/gDW for whichever
%condition is applied.
condData = readtable('biomassConditions_template.tsv', 'FileType', 'text', ...
    'Delimiter', '\t');
conditionNames = condData.Properties.VariableNames(2:end);

for c = 1:numel(conditionNames)
    conditionName = conditionNames{c};
    condModel = model; %Start each condition from the same unscaled template

    fprintf('\n--- Condition: %s ---\n', conditionName);
    for i = 1:height(condData)
        component = condData.component{i};
        target = condData.(conditionName)(i);
        %Balance the total back to 1 g/gDW on the last explicit component,
        %so carbohydrate absorbs whatever measuring the other four leaves.
        if i == height(condData)
            condModel = scaleBiomassFraction(condModel, biomassConfig, ...
                component, target, 'balanceOut', 'carbohydrate');
        else
            condModel = scaleBiomassFraction(condModel, biomassConfig, ...
                component, target);
        end
    end

    fractions = getBiomassFractions(condModel, biomassConfig);
    fprintf('carbohydrate (balanced): %.4g g/gDW\n', fractions.carbohydrate);
    fprintf('total: %.4g g/gDW\n', fractions.total);
end

%setGAM applies a known growth-associated maintenance value directly,
%rather than deriving one from data -- use fitParameters for that instead,
%against measured exchange fluxes at a fixed growth rate (see its "See
%also" cross-reference to setGAM).
before = full(model.S(strcmp(model.mets, 'atp_c'), strcmp(model.rxns, biomassConfig.biomass_rxn)));
model = setGAM(model, 35, biomassConfig.biomass_rxn, ...
    {'ATP', 'ADP', 'H2O', 'H+', 'phosphate'});
after = full(model.S(strcmp(model.mets, 'atp_c'), strcmp(model.rxns, biomassConfig.biomass_rxn)));
fprintf('\nATP coefficient in the biomass reaction: %.4g -> %.4g\n', before, after);


function model = buildBiomassTemplateModel()
%A minimal model with one pseudoreaction per biomass component, feeding
%into a top-level biomass pseudoreaction -- the layout getBiomassFractions,
%scaleBiomassFraction, scaleBiomassPseudoreaction and setGAM all expect.
model.id = 'biomassTemplateModel';
model.name = 'Template model for the biomass-composition tutorial';

model.mets = {'aa_c'; 'glycogen_c'; 'ntp_c'; 'dntp_c'; 'fa_c'; ...
    'protein_c'; 'carbohydrate_c'; 'rna_c'; 'dna_c'; 'lipid_c'; ...
    'atp_c'; 'adp_c'; 'h2o_c'; 'h_c'; 'pi_c'};
model.metNames = {'amino acid pool'; 'glycogen'; 'NTP pool'; 'dNTP pool'; ...
    'fatty acid pool'; 'protein'; 'carbohydrate'; 'RNA'; 'DNA'; ...
    'lipid_backbone'; 'ATP'; 'ADP'; 'H2O'; 'H+'; 'phosphate'};
model.metFormulas = {'C5H9NO2'; 'C6H10O5'; 'C10H12N5O10P'; 'C10H12N5O9P'; ...
    'C16H30O2'; ''; ''; ''; ''; ''; ...
    'C10H12N5O13P3'; 'C10H12N5O10P2'; 'H2O'; 'H'; 'HO4P'};
model.metCharges = [0; 0; -2; -2; 0; 0; 0; 0; 0; 0; -4; -3; 0; 1; -2];
model.comps = {'c'};
model.compNames = {'cytosol'};
model.metComps = ones(numel(model.mets), 1);

model.rxns = {'R_protein'; 'R_carb'; 'R_rna'; 'R_dna'; 'R_lipid'; 'R_biomass'};
model.rxnNames = {'protein pseudoreaction'; 'carbohydrate pseudoreaction'; ...
    'RNA pseudoreaction'; 'DNA pseudoreaction'; ...
    'lipid backbone pseudoreaction'; 'biomass pseudoreaction'};

nMets = numel(model.mets);
nRxns = numel(model.rxns);
model.S = sparse(nMets, nRxns);
met = @(id) find(strcmp(model.mets, id));
model.S(met('aa_c'), 1) = -20;         model.S(met('protein_c'), 1) = 1;
model.S(met('glycogen_c'), 2) = -5;    model.S(met('carbohydrate_c'), 2) = 1;
model.S(met('ntp_c'), 3) = -10;        model.S(met('rna_c'), 3) = 1;
model.S(met('dntp_c'), 4) = -2;        model.S(met('dna_c'), 4) = 1;
model.S(met('fa_c'), 5) = -0.08;       model.S(met('lipid_c'), 5) = 1;
model.S(met('protein_c'), 6) = -1;     model.S(met('carbohydrate_c'), 6) = -1;
model.S(met('rna_c'), 6) = -1;         model.S(met('dna_c'), 6) = -1;
model.S(met('lipid_c'), 6) = -1;
model.S(met('atp_c'), 6) = -30;        model.S(met('h2o_c'), 6) = -30;
model.S(met('adp_c'), 6) = 30;         model.S(met('h_c'), 6) = 30;
model.S(met('pi_c'), 6) = 30;

model.lb = zeros(nRxns, 1);
model.ub = repmat(1000, nRxns, 1);
model.rev = zeros(nRxns, 1);
model.c = [0; 0; 0; 0; 0; 1];
model.b = zeros(nMets, 1);
model.genes = {};
model.rxnGeneMat = sparse(nRxns, 0);
model.grRules = repmat({''}, nRxns, 1);
end


function biomassConfig = buildBiomassTemplateConfig()
biomassConfig.biomass_rxn = 'R_biomass';
biomassConfig.proton_met = 'h_c';
biomassConfig.components = { ...
    struct('name', 'protein', 'pseudoreaction_name', 'protein pseudoreaction', ...
        'mass_strategy', 'mw_minus_2h'), ...
    struct('name', 'carbohydrate', 'pseudoreaction_name', 'carbohydrate pseudoreaction', ...
        'mass_strategy', 'mw'), ...
    struct('name', 'RNA', 'pseudoreaction_name', 'RNA pseudoreaction', ...
        'mass_strategy', 'mw_minus_water'), ...
    struct('name', 'DNA', 'pseudoreaction_name', 'DNA pseudoreaction', ...
        'mass_strategy', 'mw_minus_water'), ...
    struct('name', 'lipid_backbone', 'pseudoreaction_name', 'lipid backbone pseudoreaction', ...
        'mass_strategy', 'grams')};
end
