% uih = pushbutt(h_pos,v_pos,width,height,call_st,text,inter_cbk)
% Create push button in current window
%
% h_pos		= Horizontal position (characters)
% v_pos		= Vertical position (lines), from figure top
% width  	= Width (characters).  if isempty(width),  default = 2+length(text)
% height 	= Height (characters). if isempty(height), default = 1
% call_st	= callback (string)
% text		= text
% inter_cbk		= 1-> set callback interruptible, 0-> not interruptible
%
% uih		= handle to uicontrol created
%
% version 3.3, Juan M. Rius, Jan 1997

function uih=pushbutt(h_pos,v_pos,width,height,call_st,text,inter_cbk)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;

v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top

if isempty(width), width = 2+fac*length(text); end
if isempty(height), height = 1; end

pu_pos = [h_pos*CHW v_pos*SEP width*CHW height*SEP];
uih = uicontrol(gcf,'Style','push','Position',pu_pos,'String',text,'Callback',call_st);

if nargin == 7,
	if inter_cbk == 1, set(uih,'Inter','on');
	end
end


