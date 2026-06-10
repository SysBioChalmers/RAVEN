function orangeString = printOrange(stringToPrint)
% printOrange  Print an orange-coloured string to the MATLAB Command Window.
%
% Prints stringToPrint in orange colour to the MATLAB Command Window. Only
% works if MATLAB is open with a GUI; it does not work with command-line
% MATLAB.
%
% Parameters
% ----------
% stringToPrint : char
%     string that should be printed in orange colour.
%
% Returns
% -------
% orangeString : char
%     the input string wrapped with the formatting codes for orange colour.
%
% Examples
% --------
%     printOrange(stringToPrint);

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
