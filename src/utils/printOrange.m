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
%     string that should be printed in orange colour. "\n" marks a line
%     break, as in fprintf; no other escape sequence or format directive
%     is interpreted, so text containing "%" prints as-is.
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
    orangeString = strjoin(orangeString,newline);
    %stringToPrint is arbitrary text (e.g. a model id) and may contain a
    %literal "%", which fprintf would otherwise reinterpret as a format
    %directive; the "\n" line breaks and the desktop "[\b ... ]\b"
    %highlight markers still need interpreting, so expand those explicitly
    %instead of using orangeString itself as the format string.
    orangeString = strrep(orangeString,'\n',newline);
    orangeString = strrep(orangeString,'[\b',sprintf('[\b'));
    orangeString = strrep(orangeString,']\b',sprintf(']\b'));
    fprintf('%s',orangeString);
end
end
