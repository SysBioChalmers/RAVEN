% cleanSheet
%   Cleans up an Excel sheet by removing empty rows/colums (and some other
%   checks)
%
%   raw             cell array with the data in the sheet
%   removeComments  true if commented lines (non-empty first cell in each
%                   row) should be removed (optional, default true)
%   removeOnlyCap   remove columns with captions but no other values (optional,
%                   default false)
%   removeNoCap     remove columns without captions (optional, default true)
%   removeEmptyRows remove rows with no non-empty cells (optional, default true)
%   
%   raw             cleaned version
%   keptRows        indexes of the kept rows in the original structure
%   keptCols        indexes of the kept columns in the original structure
%
% Usage: [raw,keptRows,keptCols]=cleanSheet(raw,removeComments,removeOnlyCap,...
%               removeNoCap,removeEmptyRows)

function [raw,keptRows,keptCols]=cleanSheet(raw,removeComments,removeOnlyCap,removeNoCap,removeEmptyRows)
if nargin<2
    removeComments=true;
end
if nargin<3
    removeOnlyCap=false;
end
if nargin<4
    removeNoCap=true;
end
if nargin<5
    removeEmptyRows=true;
end

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
