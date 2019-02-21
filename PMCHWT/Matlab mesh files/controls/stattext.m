% uih = stattext(h_pos,v_pos,text,fg_color,bg_color,width)
% Create static text in current window
%
% h_pos		= Horizontal position (characters)
% v_pos		= Vertical position (lines), from figure top
% text		= text string
% fg_color	= foreground color [R G B], if fg_color=[] or nargin<4, default = [0 0 0]
% bg_color	= background color [R G B], if bg_color=[] or nargin<5, default = get(gcf,'Color')
%
% uih	= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Sept 1996

function uih=stattext(h_pos,v_pos,text,fg_color,bg_color)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top

if nargin<5, bg_color = color;
elseif isempty(bg_color), bg_color = color;
end

if nargin<4, fg_color = [0 0 0];
elseif isempty(fg_color), fg_color = [0 0 0];
end

tit_pos = [h_pos*CHW v_pos*SEP (2+fac*length(text))*CHW CHH];
uih = uicontrol(gcf,'Style','Text','String',text,'Position',tit_pos,...
	'BackGroundColor',bg_color,'ForeGroundColor',fg_color,'HorizontalAlignment','Left');


