function tbl=readKEGGTable(gzFile)
% readKEGGTable  Read a gzipped KEGG relational table (published by
% raven-toolbox as gzipped TSV, see raven-toolbox's
% docs/maintenance/kegg_data_format.md) into a table with every column
% forced to char cellstr (matching RAVEN's own cell-array-of-char
% convention, regardless of the caller's default text-import type).
%
% Parameters
% ----------
% gzFile : char
%     path to the gzipped TSV file, e.g. '.../kegg118_ko_reaction.tsv.gz'.
%
% Returns
% -------
% tbl : table
%     the table, with every column as a cellstr.

plainFile=gzFile(1:end-3); %strip the trailing '.gz'
if ~isfile(plainFile)
    gunzip(gzFile);
end
opts=detectImportOptions(plainFile,'FileType','text','Delimiter',char(9));
opts=setvartype(opts,opts.VariableNames,'char');
tbl=readtable(plainFile,opts);
end
