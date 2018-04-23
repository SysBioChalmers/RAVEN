function strings=cleanBadChars(strings)
% cleanBadChars
%   Converts characters which are illegal in SBML to their entity reference
%
%   strings     a string or cell array
%
%   %strings    cleaned string or cell array
%
%   The characters <>&'" are illegal in SBML (in most cases). They are here
%   converted to their entity references. Empty or non-string cells are
%   ignored.
%
%   Usage: strings=cleanBadChars(strings)
%
%   Rasmus Agren, 2013-08-03
%

%For simplicity
wasChar=false;
if ischar(strings)
    strings={strings};
    wasChar=true;
end

%Get the ids which are cell strings
I=cellfun(@isstr,strings);

%Note that the order here matters
strings(I)=strrep(strings(I),'&','&amp;');
strings(I)=strrep(strings(I),'<','&lt;');
strings(I)=strrep(strings(I),'>','&gt;');
strings(I)=strrep(strings(I),'"','&quot;');
strings(I)=strrep(strings(I),'''','&apos;');
strings(I)=strrep(strings(I),'%','&#37;');

if numel(strings)==1 && wasChar==true
    strings=strings{1};
end
end
