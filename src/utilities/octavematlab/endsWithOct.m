function TF = endsWithOct(STR,PAT)
    % TODO: doc
    % TODO: implement ignore case parameter
if ~isOctave && verLessThan('matlab','9.1') % function introduced R2016b
    TF = contains(STR,PAT)
else
    specialchars = '[]{}().,*?!+=@$|\';
    specialchars = num2cell(specialchars);
    specialchars = ['(\' strjoin(specialchars,')|(\') ')'];
    PAT = regexprep(PAT,specialchars,'\\$0');
    PAT = [PAT '$'];
    if ischar(PAT)
        PAT={PAT};
    end
    if ischar(STR)
        STR={STR};
    end
    TF = zeros(numel(STR),1);
    for i = 1:numel(PAT)
        TF(~cellfun(@isempty,regexp(STR,PAT{i}))) = 1;
    end
    TF = logical(TF);
%    if (numel(TF) == 1)
%        TF = TF {1};
%    end
end
end
