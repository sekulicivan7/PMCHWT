% uih = popup(h_pos,v_pos,var,value,text,tx_size,ed_size)
% Create popup box in current window
% Changes global variable or gcf's UserData
%
% h_pos	= Horizontal position (characters)
% v_pos	= Vertical position (lines), from figure top
% var	= Global variable name (string), or gcf's UserData index (number)
%	  Strings are stored in gcf's UserData as handle to invisible StaticText uicontrol
%	  Data must be recovered as get(gud(gcf,var),'String');
% value	= matrix of strings, containing menu options
% text	= text string
% tx_size = (optional) size of text box (characters)
% ed_size = (optional) size of edit box (characters)
%
% uih	= handle to uicontrol created
%
% version 3.2, Juan M. Rius, Oct 1996

function uih=popup(h_pos,v_pos,var,value,text,tx_size,ed_size)

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;

v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top
        
if nargin < 7, ed_size = 5*(type+1);		% Default size edit box:
elseif isempty(ed_size), ed_size = 5*(type+1);	% numeric->5, string->10
end

if nargin < 6, tx_size = 20;			% Default size text box = 20
elseif isempty(ed_size), tx_size = 20;
end
   
ed_pos = [(h_pos+tx_size+1)*CHW SEP*v_pos ed_size*CHW CHH];
tx_pos = [h_pos*CHW SEP*v_pos tx_size*CHW CHH];

[M,N]=size(value);

% Create string with items separated by '|'
string = value(1,:); for m=2:M, string = [string '|' deblank(value(m,:))]; end

if isstr(var),
	call_st=[var '=deblank(gud(gco,get(gco,''Value''),0));'];
else,
	sud(gcf,uicontrol('Style','Text','String',value(1,:),'Visible','off'),var);
	call_st=['set(gud(gcf,' num2str(var) '),''String'',deblank(gud(gco,get(gco,''Value''),0)));'];
end;

ui_ed = uicontrol(gcf,'Style','popup','Position',ed_pos,'HorizontalAlignment','Right',...
			'BackGroundColor',[0.5 0.5 0.5],'ForeGroundColor',[1 1 1],...
			'UserData',value,'String',string,'Callback',call_st);

ui_tx = uicontrol(gcf,'Style','Text','String',text,'Position',tx_pos,...
			'BackGroundColor',color,'ForeGroundColor',[0 0 0],...
			'HorizontalAlignment','Right');

uih = [ui_ed ui_tx];

