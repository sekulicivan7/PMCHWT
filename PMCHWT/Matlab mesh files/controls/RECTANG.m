% rec = rectang(h_pos,v_pos,width,height,color)
%
% Draws rectangle in new axis, in "character" units 
%
% h_pos  = horizontal position in characters
% v_pos  = vertical position in characters (from top)
% width  = width in characters
% height = height in characters
% color  = rectangle color
%
% rec = handle to patch object
%
% Juan M. Rius
% July 1996

function rec = rectang(h_pos,v_pos,width,height,color)

[CHH,CHW,SEP,back,H_SIZE,V_SIZE]=init_win;
v_pos = V_SIZE-v_pos;
axe_e = axes('Position',[0 0 0.01 0.01]);

set(gca,'Units','pixels','Position',[h_pos*CHW v_pos*SEP width*CHW height*SEP],...
	'Box','on','Xticklabels',[],'Yticklabels',[],'drawmode','fast',...
	'XColor','k','YColor','k','Xtick',[0 1],'Ytick',[0 1]);

rec = patch([0 1 1 0],[0 0 1 1],color,'EdgeColor',color,'EraseMode','background');

pix_x = 1/(width*CHW);
pix_y = 2/(height*SEP);
axis([-pix_x 1+pix_x -pix_y  1+pix_y]);

