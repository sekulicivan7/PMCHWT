% uih = frame(h_pos,v_pos,width,height,fg_color,bg_color)
% Create frame in current window
%
% h_pos		= Horizontal position (characters)
% v_pos		= Vertical position (lines), from figure top
% width 	= Width (characters)
% height	= Height (characters)
% fg_color	= foreground color [R G B], if fg_color=[] or nargin<5, default = [0 0 0]
% bg_color	= background color [R G B], if bg_color=[] or nargin<6, default = get(gcf,'Color')
%
% uih		= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Sept 1996

function uih=frame(h_pos,v_pos,width,height,fg_color,bg_color)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top        

if nargin<6, bg_color = color;
elseif isempty(bg_color), bg_color = color;
end

if nargin<5, fg_color = [0 0 0];
elseif isempty(fg_color), fg_color = [0 0 0];
end

fr_pos = [h_pos*CHW v_pos*SEP width*CHW height*SEP];
uih = uicontrol(gcf,'Style','frame','Position',fr_pos,'BackGroundColor',bg_color,'ForeGroundColor',fg_color);


