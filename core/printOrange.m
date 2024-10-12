function orangeString = printOrange(stringToPrint)
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
    orangeString = ['[\b' stringToPrint,']\b'];
else
    orangeString = stringToPrint;
end
if nargout < 1
    % Wrap text to command window size
    sz = get(0, 'CommandWindowSize');
    orangeString = textwrap({orangeString},sz(1));
    orangeString = strjoin(orangeString,'\n');
    fprintf(orangeString);
end
end
