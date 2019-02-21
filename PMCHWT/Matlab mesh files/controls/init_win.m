% [CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win(win_h,H_SIZE,V_SIZE,title,color,type)
%
% With only input arguments:
% 	Initialize figure for uicontrols
% 	The default is character height 20 pixels for 1024x768 display, proportional for others
% 	This default works well with 9 point Arial uicontrol font
%
%	win_h  = Handle to window
% 	H_SIZE = Horizontal size of figure (characters)
% 	V_SIZE = Vertical size of figure (lines)
% 	title  = Window name (string) (optional). Default = 'Input parameters'
% 	color  = Background color (optional). Default = [0.76,0.76,0.76].
%	type   = 'window' or 'plot' (optional). Default 'window'.
%
% With only output arguments
% 	Get figure information necessary for uicontrol functions
%
% 	CHH  = Character height (pixels)
% 	CHW  = Character width (pixels)
%	fac  = Reduction factor in character width due to narrow font
% 	SEP  = Interline separation (pixels)
% 	color = Background color of figure
% 	H_SIZE = Horizontal size of figure (characters)
% 	V_SIZE = Vertical size of figure (lines)
%
% version 3.3,	Juan M. Rius, Jan 1997
% Version 3.3: 	a) can be called from a callback function (interruptible)
%		b) deletes the figure
% BE CAREFUL!! If init_win is called from callback, you must set 'visible off'
% on creation of the figure. Otherwise set(fig,'vis','off') will change current figure
% and uicontrols will be drawn on the wrong figure

function [CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win(win_h,H_SIZE,V_SIZE,title,color,type)

screen = get(0, 'ScreenSize');
CHH = 20*(screen(3)/1024); CHW=CHH/2;	% Character size = 20 pixels height for 1024x768
SEP = CHH*1.25;				% Interline separation (pixels)
fac = 0.8;

if nargin == 0,
	win_h = gcf;
	wpos=get(win_h,'Position');
	H_SIZE=wpos(3)/CHW; V_SIZE=wpos(4)/SEP;
	color = get(win_h,'Color');
elseif nargout==0,

	delete(get(win_h,'Children'));

	if nargin < 6, type = 'window';
	end

	if nargin < 5, color = [0.76,0.76,0.76];
	end;

	if nargin < 4, title = 'Input parameters:';
	end

	BORDE = 5;	% Thickness of MATLAB window border

	wpos=get(win_h,'Position');
   wpos(1) = screen(3)-BORDE-H_SIZE*CHW-40; wpos(2) = screen(4)-BORDE-(V_SIZE+1)*SEP;
	wpos(3) = H_SIZE*CHW;                 wpos(4) = V_SIZE*SEP;

	if strcmp(type,'window'),
      set(win_h,'Position',wpos,'Resize','off','MenuBar','none',...
			'NextPlot','new',...
         'Color',color,'Numbertitle','Off','Name',title);
	elseif strcmp(type,'plot'),
			set(win_h,'Position',wpos,'Color',color,'Numbertitle','Off','Name',title);
	else error('type argument must be window or plot');
	end

else error('init_win must be called with either 0 input parameters or 0 output parameters');
end
