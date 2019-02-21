% butt = dlgwin(type,title,text,pushst1,pushst2,pushst3,pushst4,pushst5,pushst6,pushst7)
% Displays a dialog window with text and pushbuttons,
% waits for click in a pushbutton and returns number of pushbutton clicked
% dlgwin() can be called from a callback function (interruptible)
%
% Output:
% butt = Number of pushbutton clicked
%
% Input:
% type  = 1-> Question, wait for answer; 0-> Information, continue processing
%	  Type 0 (information) do not return number of pushbutton clicked
%	  Type 1 (question) requieres INTERRUPTIBLE calbback calling dlgwin
% title	= title of window
% text	= text in window (one line)
% pushst1 = string in 1st pushbutton
% pushst2 = string in 2nd pushbutton
% ... etc
%
% version 3.3, Juan M. Rius, Jan 1997

function butt = dlgwin(type,title,text,pushst1,pushst2,pushst3,pushst4,pushst5,pushst6,pushst7)

DEF_WIDTH = 20;		% Initial width, will be changed latter

% BE CAREFUL!! If init_win is called from callback, you must set 'visible off'
% on creation of the figure. Otherwise set(fig,'vis','off') will change current figure
% and uicontrols will be drawn on the wrong figure

fig = figure('Vis','off');

init_win(fig,DEF_WIDTH,3,title);

% Count width of text and pushbuttons
widtx = 0.9*length(text);
widbut = 2;
for i=1:nargin-3,
	st = eval(sprintf('pushst%d',i));
	widbut = widbut + length(st) + 2;
end

% Compute position of text and first pushbutton
sh = (widtx - widbut)/2;

if sh>=0,	stattext(1,1,text);
		x = sh;		% Position of pushbutton
else,		stattext(-sh,1,text);
		x = 1;		% Position of pushbutton
end

pbs = [];			% List oh handles to pushbuttons
for i=1:nargin-3,
	st = eval(sprintf('pushst%d',i));
	if type, cbk = ';';	% No callback is necessary
	else cbk = 'close(gcf)';% Information: close figure
	end
	p = pushbutt(x,2.5,[],[],cbk,st);
	pbs = [pbs p];		% Save handle to detect click (later)
	set(p,'UserData',i);	% Number of pushbutton (to return)
	x = x + length(st) + 2; % Update position for next pushbutton
end

% Center window and modify width according to length of text and pushbuttons
screen = get(0, 'ScreenSize');
wpos = get(fig,'Position');
wpos(3) = wpos(3) * max(widbut,widtx)/DEF_WIDTH;
wpos(2) = (screen(3) - wpos(3))/2;
wpos(1) = (screen(4) - wpos(4))/2;
set(fig,'Position',wpos);

% Display figure and wait for click
set(fig,'Vis','on');

if type,
	while ~any(get(fig,'CurrentObject')==pbs),
		waitforbuttonpress;
	end

	% Return number of pushbutton
	butt = get(get(fig,'CurrentO'),'UserData');

	% Uncomment this line to return the string in the pushbutton
	% butt = get(get(fig,'CurrentO'),'Str');
	delete(fig);
end

