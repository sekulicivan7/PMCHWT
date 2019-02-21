% addcallb(uih,cmd_before,cmd_after)
% Add commands to callback of uicontrol
%
% uih		= handle to uicontrol
% cmd_before	= command string to be executed before
% cmd_after	= command string to be executed after
%
% version 3.0, Juan M. Rius, Oct 1996

function addcallb(uih,cmd_before,cmd_after)

if isstr(cmd_before) & isstr(cmd_after),
	set(uih,'Callback',[cmd_before ';' get(uih,'Callback') cmd_after ';']);
else	error('command argument is not string');
end
