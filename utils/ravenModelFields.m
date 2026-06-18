function reg = ravenModelFields()
% ravenModelFields  Declarative registry of all RAVEN model struct fields.
%
% Returns a struct array. Each element describes one field:
%   .name    - field name (char)
%   .type    - entity type: 'rxn', 'met', 'gene', or 'comp' (char)
%   .default - default value when padding new entries
%
% Only 1-D (vector-per-entity) fields are listed here. Fields that require
% custom 2-D indexing or value-remapping are excluded and must be handled
% explicitly by callers:
%
%   S          - sparse (mets x rxns); use S(:,perm) for rxns, S(perm,:) for mets
%   rxnGeneMat - sparse (rxns x genes); use rxnGeneMat(perm,:) or rxnGeneMat(:,perm)
%   b          - double (mets x N), may be 2-D; use b(perm,:)
%   metComps, rxnComps, geneComps - when permuting comps, their VALUES must
%       be remapped (not the arrays themselves permuted)

rows = {
% name                   type      default
  'rxns',                'rxn',    '';
  'lb',                  'rxn',    0;
  'ub',                  'rxn',    1000;
  'rev',                 'rxn',    false;
  'c',                   'rxn',    0;
  'grRules',             'rxn',    '';
  'rxnNames',            'rxn',    '';
  'eccodes',             'rxn',    '';
  'subSystems',          'rxn',    '';
  'rxnFrom',             'rxn',    '';
  'rxnNotes',            'rxn',    '';
  'rxnReferences',       'rxn',    '';
  'pwys',                'rxn',    '';
  'equations',           'rxn',    '';
  'rxnComps',            'rxn',    NaN;
  'rxnConfidenceScores', 'rxn',    NaN;
  'rxnDeltaG',           'rxn',    NaN;
  'rxnScores',           'rxn',    NaN;
  'spontaneous',         'rxn',    false;
  'rxnMiriams',          'rxn',    {};
  'mets',                'met',    '';
  'metComps',            'met',    0;
  'metNames',            'met',    '';
  'metFormulas',         'met',    '';
  'inchis',              'met',    '';
  'metSmiles',           'met',    '';
  'metNotes',            'met',    '';
  'metFrom',             'met',    '';
  'metCharges',          'met',    NaN;
  'metDeltaG',           'met',    NaN;
  'unconstrained',       'met',    0;
  'metMiriams',          'met',    {};
  'genes',               'gene',   '';
  'geneShortNames',      'gene',   '';
  'geneFrom',            'gene',   '';
  'proteins',            'gene',   '';
  'geneComps',           'gene',   NaN;
  'geneMiriams',         'gene',   {};
  'comps',               'comp',   '';
  'compNames',           'comp',   '';
  'compOutside',         'comp',   '';
  'compMiriams',         'comp',   {};
};

n = size(rows, 1);
for i = 1:n
    reg(i).name    = rows{i,1};
    reg(i).type    = rows{i,2};
    reg(i).default = rows{i,3};
end
end
