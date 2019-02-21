% uih = inputbox(h_pos,v_pos,var,value,text,tx_size,ed_size)
% Create input box in current window
% Changes global variable or gcf's UserData
%
% h_pos	= Horizontal position (characters)
% v_pos	= Vertical position (lines), from figure top
% var	= variable name (string)
% var	= Global variable name (string), or gcf's UserData index (number)
%	  Strings are stored in gcf's UserData as handle to invisible StaticText uicontrol
%	  Data must be recovered as get(gud(gcf,var),'String');
% value	= initial contents (string or number). May be void: [] or ''
% text	= text string
% tx_size = (optional) size of text box (characters)
% ed_size = (optional) size of edit box (characters)
%
% uih	= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Oct 1996

function uih=inputbox(h_pos,v_pos,var,value,text,tx_size,ed_size)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;

type = isstr(value);			% type = 1 if is string
if ~type,
	[M,N] = size(value);
	if N*M>1,			% value is a matrix
		tmp = value;
		value = '[';
		for m=1:M,
			for n=1:N,  value = [value num2str(tmp(m,n)) ' ']; end
			value = [value(1:length(value)-1) '; '];
		end
		value = [value(1:length(value)-2) ']'];
	else	value = num2str(value); % value is a number
	end
end

v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top

if nargin < 7, ed_size = 5*(type+1);		% Default size edit box:
elseif isempty(ed_size), ed_size = 5*(type+1);	% numeric->5, string->10
end

if nargin < 6, tx_size = 20;			% Default size text box = 20
elseif isempty(ed_size), tx_size = 20;
end
   
ed_pos = [(h_pos+tx_size+1)*CHW SEP*v_pos 1.1*ed_size*CHW CHH];
tx_pos = [h_pos*CHW SEP*v_pos tx_size*CHW CHH];

if isstr(var),
	if type==1,call_st = [ var '=get(gco,''String'');'];
	else	   call_st = [ var '=str2num(get(gco,''String''));'];
	end;
else,
	if type==1,sud(gcf,uicontrol('Style','Text','String',value,'Visible','off'),var);
		   call_st = ['set(gud(gcf,' num2str(var) '),''String'',get(gco,''String''));' ];
	else	   call_st = ['sud(gcf,str2num(get(gco,''String'')),' num2str(var) ');'];
	end;
end;

ui_ed = uicontrol(gcf,'Style','Edit','Position',ed_pos,'HorizontalAlignment','Left',...
			'BackGroundColor',[0.5 0.5 0.5],'ForeGroundColor',[1 1 1],...
			'String',value,'Callback',call_st);

ui_tx = uicontrol(gcf,'Style','Text','String',text,'Position',tx_pos,...
			'BackGroundColor',color,'ForeGroundColor',[0 0 0],...
			'HorizontalAlignment','Right');

uih = [ui_ed ui_tx];

