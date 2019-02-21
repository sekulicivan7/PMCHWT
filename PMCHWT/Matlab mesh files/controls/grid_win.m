% grid_win
% Draw lines-characters grid in current window
%
% Version 3.2, 	October 1996, J.M.Rius

function grid_win

[CHH,CHW,fac,SEP,color,H_SIZE,V_SIZE]=init_win;

axes('position',[0 0 1 1])
set(gca,'Xlim',[0 H_SIZE])
set(gca,'Ylim',[0 V_SIZE])
set(gca,'Xtick',0:50)
set(gca,'Ytick',0:20)
set(gca,'GridLineStyle','-')
if any(color), set(gca,'xcolor','k'); set(gca,'ycolor','k'); end
grid on

for i=1:V_SIZE,
	text(0.25,V_SIZE+0.5-i,num2str(i),'color','r','fontsize',8,'fontname','helvetica');
	text(H_SIZE-1.5,V_SIZE+0.5-i,num2str(i),'color','r','fontsize',8,'fontname','helvetica');
end

for i=1:H_SIZE,
	text(i-0.5,0.25,num2str(i),'color','r','fontsize',6,'fontname','helvetica','rotation',90);
	text(i-0.5,V_SIZE-0.7,num2str(i),'color','r','fontsize',6,'fontname','helvetica','rotation',90);
end
