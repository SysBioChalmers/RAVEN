function printOrange(stringToPrint)
% printOrange
%   Print orange-colored stringToPrint to the MATLAB Command Window. Only
%   if MATLAB is open with GUI, does not work with command-line MATLAB.
%
% Input: 
%   stringToPrint   string that should be printed in orange color
%
% Usage: printOrange(stringToPrint)
try useDesktop = usejava('desktop'); catch, useDesktop = false; end
if useDesktop
    fprintf(['[\b' stringToPrint,']\b'])
else
    fprintf(stringToPrint)
end
end
