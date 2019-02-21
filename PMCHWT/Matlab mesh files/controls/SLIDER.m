% uih = slider(h_pos,v_pos,width,height,var,value,title)
% Create horizontal slider in current window
% String with current value is editable!!
% Changes global variable or gcf's UserData
% 
% h_pos	 = Horizontal position (characters)
% v_pos	 = Vertical position (lines), from figure top
% width  = Width (characters).  if isempty(width),  default = 25
% height = Height (characters). if isempty(height), default = 1.25
% var	 = Global variable name (string), or gcf's UserData index (number)
% value	 = [initial_value, min_value, max_value], initial, maximum and minimum values
% title	 = Title (string)
%
% uih = [uih_slider uih_min uih_max uih_val uih_tit] handles of uicontrols created
%
% version 3.2, Juan M. Rius, Oct 1996

function uih=slider(h_pos,v_pos,width,height,var,value,title)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top

if isempty(width),  width = 25; end
if isempty(height), height = 1.25; end

ed_pos = [h_pos*CHW v_pos*SEP width*CHW height*SEP-CHH];
smin = num2str(value(2));	lmin = (fac*length(smin)+2)*CHW;
smax = num2str(value(3));	lmax = (fac*length(smax)+2)*CHW;
tit  = num2str(value(1));	ltex = (fac*length(tit) +2)*CHW;
				lnom = (fac*length(title)+2)*CHW;

uih_min = uicontrol(gcf,'Style','Text','BackgroundColor',color,...
	'Position',[ed_pos(1)                ed_pos(2) lmin CHH],'String',smin);
	
uih_max = uicontrol(gcf,'Style','Text','BackgroundColor',color,...
	'Position',[ed_pos(1)+ed_pos(3)-lmax ed_pos(2) lmax CHH],'String',smax);

if isstr(var),
	call_ed=[ var '=str2num(get(gco,''String'')); ',...
		 'set(gud(gco),''Value'',' var '); ',...
		];
else,
	call_ed=['sud(gcf,str2num(get(gco,''String'')),' num2str(var) '); ' ...
		 'set(gud(gco),''Value'',gud(gcf,' num2str(var) ')); ',...
		];
end;

uih_val = uicontrol(gcf,'Style','Edit','BackgroundColor',color,...
	'CallBack',call_ed,'String',num2str(value(1)),...
	'Position',[ed_pos(1)+lnom+(ed_pos(3)-lnom-ltex)/2 ed_pos(2)+height*SEP-CHH ltex CHH],'String',tit);

uih_tit = uicontrol(gcf,'Style','Text','BackgroundColor',color,...
	'Position',[ed_pos(1)+(ed_pos(3)-lnom-ltex)/2 ed_pos(2)+height*SEP-CHH lnom CHH],'String',title);
	
ed_pos(1) = ed_pos(1) + lmin; ed_pos(3) = ed_pos(3) - lmin - lmax;

if isstr(var),
	call_st=[var '= get(gco,''Value'');'...
		 'set(gud(gco),''String'',num2str(' var '));'];
else,
	call_st=['sud(gcf,get(gco,''Value''),' num2str(var) '); ' ...
		 'set(gud(gco),''String'',num2str(gud(gcf,' num2str(var) '))); '...
		];
end

uih_slider = uicontrol(gcf,'Style','Slider','Position',ed_pos,...
	'Min',value(2),'Max',value(3),'Value',value(1),'UserData',uih_val,...
	'CallBack',call_st);
sud(uih_val,uih_slider);
uih = [uih_slider uih_val uih_min uih_max uih_tit];

