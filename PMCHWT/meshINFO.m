

geom	= 'prism_Ivan';	

     	param = [0.1 7 0.1 7 0.1 7];
	   
%%*************************************************************************************************************
if (length(param)==1),
    %%% gid object file
    p_obj.gid.file = geom;
    m_gidmesh_Ivan;
else,
    %%% parameters
    obj= feval(geom,param);
    obj = get_edge_Ivan(obj);
end;

topol2(1,:)=obj.topol(1,:);
topol2(2,:)=obj.topol(3,:);    %%changing the order of vertices in triangle to be counterclockwise!!!!
topol2(3,:)=obj.topol(2,:);

trian2(1,:)=obj.trian(1,:);
trian2(2,:)=obj.trian(3,:);
trian2(3,:)=obj.trian(2,:);

vertex=obj.vertex.';
topol2=topol2.';
trian2=trian2.';

  dlmwrite('coord.txt',vertex,'delimiter','\t','precision',10)
  dlmwrite('topol.txt',topol2,'delimiter','\t')
  dlmwrite('trian.txt',trian2,'delimiter','\t')
  
clear

