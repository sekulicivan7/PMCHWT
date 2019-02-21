% uih = checkbox(h_pos,v_pos,var,value,text,control_dep)
% Create checkbox in current window
% Changes global variable or gcf's UserData
%
% h_pos		= Horizontal position (characters)
% v_pos		= Vertical position (lines), from figure top
% var	= Global variable name (string), or gcf's UserData index (number)
% value		= initial value (1->checked, 0->unchecked)
% text		= text
% control_dep 	= Optional argument: controls with visibility dependent from this box
%		[handle_1 flag_1 handle_2 flag_2 handle_3 flag_3 ...]
%		if flag_i==1, handle_i is visible when this box is   checked
%		if flag_i==0, handle_i is visible when this box is unchecked
%
% uih	= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Oct 1996

function uih=checkbox(h_pos,v_pos,var,value,text,control_dep)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top

if nargin<6, control_dep = [];
end

if isstr(var),
	call_st=[ var '=' 'get(gco,''Value'')==1; ',...
		 'visibil(gco); '...
		];
else,
	call_st=[ 'sud(gcf,get(gco,''Value'')==1,' num2str(var) '); ',...
		  'visibil(gco); '...
		];
end;

cb_pos = [h_pos*CHW v_pos*SEP (fac*length(text)+2)*CHW SEP];

uih = uicontrol(gcf,'Style','checkbox','Position',cb_pos,'String',text,'UserData',control_dep,...
	'BackGroundColor',color,'ForeGroundColor',[0 0 0],'Value',value,'Callback',call_st);

visibil(uih);


