% cleanSheet  Clean up an Excel sheet.
%
% Removes empty rows/columns (and performs some other checks).
%
% Parameters
% ----------
% raw : cell
%     cell array with the data in the sheet.
% removeComments : logical, optional
%     true if commented lines (non-empty first cell in each row) should be
%     removed (default true).
% removeOnlyCap : logical, optional
%     remove columns with captions but no other values (default false).
% removeNoCap : logical, optional
%     remove columns without captions (default true).
% removeEmptyRows : logical, optional
%     remove rows with no non-empty cells (default true).
%
% Returns
% -------
% raw : cell
%     cleaned version.
% keptRows : double
%     indices of the kept rows in the original structure.
% keptCols : double
%     indices of the kept columns in the original structure.
%
% Examples
% --------
%     [raw, keptRows, keptCols] = cleanSheet(raw, removeComments, ...
%         removeOnlyCap, removeNoCap, removeEmptyRows);
function [raw,keptRows,keptCols]=cleanSheet(raw,varargin)
p=parseRAVENargs(varargin, {'removeComments',true; 'removeOnlyCap',false; 'removeNoCap',true; 'removeEmptyRows',true});
removeComments=p.removeComments; removeOnlyCap=p.removeOnlyCap; removeNoCap=p.removeNoCap; removeEmptyRows=p.removeEmptyRows;

keptRows=1:size(raw,1);
keptRows=keptRows(:);
keptCols=1:size(raw,2);
keptCols=keptCols(:);

%First check that it's a cell array. If a sheet is completely empty, then
%raw=NaN
if iscell(raw)
    %Clear cells which contain only white spaces or NaN. This could happen
    %if you accidentally inserted a space for example. I don't know how NaN
    %could occur after switching to Apache POI, but I clear it to be sure
    whites=cellfun(@wrapperWS,raw);
    raw(whites)={[]};
    
    %Find the rows that are not commented. This corresponds to the first
    %row and the ones which are empty in the first column
    if removeComments==true
        keepers=cellfun(@isempty,raw(:,1));
        keepers(1)=true;
        raw=raw(keepers,:);
        keptRows=keptRows(keepers);
    end
    
    %Remove columns that don't have string headers. If you cut and paste a
    %lot in the sheet there tends to be columns that are empty
    if removeNoCap==true
        I=cellfun(@isstr,raw(1,:));
        raw=raw(:,I);
        keptCols=keptCols(I);
    end
    
    %Remove columns which are empty except for header
    if removeOnlyCap==true
        I=~all(cellfun(@isempty,raw(2:end,:)),1);
        raw=raw(:,I);
        keptCols=keptCols(I);
    end
    
    %Check if there are any rows that are all empty. This could happen if
    %it reads too far or if the user has inserted them for layout reasons.
    %Remove any such rows
    if removeEmptyRows==true
        I=~all(cellfun(@isempty,raw),2);
        raw=raw(I,:);
        keptRows=keptRows(I);
    end
else
    raw={[]};
end

%Checks if something is all white spaces or NaN
    function I=wrapperWS(A)
        if isnan(A)
            I=true;
        else
            %isstrprop gives an error if boolean
            if islogical(A)
                I=false;
            else
                %I have to make this check, since otherwise it will pick up
                %on doubles for which the ASCII representation is a white
                %space character
                if isnumeric(A)
                    I=false;
                else
                    I=all(isstrprop(A,'wspace'));
                end
            end
        end
    end
end
