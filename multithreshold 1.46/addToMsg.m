% ---------------------------------------------------------- addToMsg
function addToMsg(message,append, warning)
% 'message' is posted to the message board
% if append==1, message is appended to the current message
global experiment expGUIhandles

if nargin<3
    % this is not a warning screen
    warning=0;
end

if append
msg=get(expGUIhandles.textMSG,'string');
[r c]=size(msg);
    if length(message)<=c
        y=[message blanks(c-length(message))];
        msg(r+1,:)=y;
    else
        msg=message;
    end
else
    msg=message;
end

try
    set(expGUIhandles.textMSG,'string', msg,'fontSize', experiment.msgFontSize)
    if warning
        % flash red to signal a warning
        set(expGUIhandles.textMSG,'backgroundcolor','r', 'ForegroundColor', 'w'	)
    else
        set(expGUIhandles.textMSG,'backgroundcolor','w', 'ForegroundColor', 'b')
    end
catch
    error(message)
end