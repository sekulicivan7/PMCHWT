% uih = radiobut(h_pos,v_pos,var,value,text,control_dep)
% Create array of radio buttons in current window
% Changes global variable or gcf's UserData
%
% h_pos	= Array of horizontal positions (characters)
% v_pos	= Array of vertical positions (lines), from figure top
% var	= Global variable name (string), or gcf's UserData index (number)
% value	= [initial_value, value_1, value_2, value_3, ...]
% text	= Array of text strings
% control_dep = Optional argument: controls with visibility dependent from this buttons
%		[handle_1 value_1_1 value_1_2 ... handle_2 value_2_1 value_2_2 ... ]
%		value_i_* are the possible values of 'var' for which handle_i is visible
%		for the other values, handle_i is invisible
%
% uih	= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Oct 1996

function uih=radiobut(h_pos,v_pos,var,value,text,control_dep)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top        

if nargin<6, control_dep = [];
end

N = length(h_pos);	% Number of buttons in the array
mat = [zeros(1,N+1) control_dep];
mat(1) = N;		% mat = [N h1 h2 h3 ... hN, control_dep];

for i=1:N,
	tx = deblank(text(i,:));
	cb_pos = [h_pos(i)*CHW v_pos(i)*SEP (fac*length(tx)+2)*CHW SEP];
	mat(i+1) = uicontrol(gcf,'Style','radiobutton','Position',cb_pos,...
	'String',tx,'BackGroundColor',color,'ForeGroundColor',[0 0 0],...
	'Max',value(i+1),'Min',-1,'Value',-1);
end

for i=1:N,
	if isstr(var),
		call_st = ['setradio(gco); ',...
			   var '= get(gco,''Value''); ',...
			   'visibil(gco); ',...
			  ];
	else,
		call_st = ['setradio(gco); ' ...
			   'sud(gcf,get(gco,''Value''),' num2str(var) '); ' ...
			   'visibil(gco); ' ...
			  ];
	end
	% Store: call_st, UserData, value. Cannot be done before because UserData contains button handles
	if value(1) == value(i+1), set(mat(i+1),'CallBack',call_st,'UserData',mat,'Value',get(mat(i+1),'Max'));
	else 			   set(mat(i+1),'CallBack',call_st,'UserData',mat,'Value',get(mat(i+1),'Min'));
	end
end

uih = mat(2:N+1);

visibil(uih);

