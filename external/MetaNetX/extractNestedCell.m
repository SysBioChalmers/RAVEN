function [refID,extractedID] = extractNestedCell(referenceCell,nestedCell)
%editMNXrefFields  Extract values in nested cells with the values in the reference cell
%
%   referenceCell        cell array, in which each cell is linked to the
%                        corresponding individual cells in nestedCell 
%               
%   nestedCell           cell array which contains the nested cells in some
%                        or all of the cells
%               
%   E.g. A is a cell array containing cells with unique metMetaNetXID, while
%        B is a cell array containing nested 1x2 cells with metChEBIID
%        Each cell in B is matched to the cell in A i.e. ID in A(1,1) is
%        linked to IDs in nested cell B(1,1)
%        A will be the referenceCell, while B will be the nestedCell
%
%   Note that RAVEN function flattenCell.m is used
%
% Usage: [refID,extractedID] = extractNestedCell(referenceCell,nestedCell)
%
% Cheng Wei Quan (Eiden), 2020-05-14

temp(:,1) = referenceCell;
temp(:,2) = nestedCell;
empties = find(cellfun('isempty',temp(:,2)));
temp(empties,:) = [];
extractNested = flattenCell(temp(:,2));
totalIDcount = find(~cellfun('isempty',extractNested));
refID{size(totalIDcount,1),1} = [];
extractedID{size(totalIDcount,1),1} = [];

for i = 1:size(extractNested,1)
    count = find(~cellfun('isempty',extractNested(i,:)));
    for j = 1:size(count,2)
        arrayidx = find(cellfun('isempty',refID),1);
        refID(arrayidx) = temp(i,1);
        arrayidx2 = find(cellfun('isempty',extractedID),1);
        extractedID(arrayidx2) = extractNested(i,j);
    end
end

[extractedID,idx] = sort(extractedID);
refID = refID(idx);

end
