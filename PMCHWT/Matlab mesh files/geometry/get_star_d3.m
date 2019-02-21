% Procedure in order to get all the variables related to STAR basis
% functions
%
% There are (Nt-1)xSTAR basis functions
% Unless differently stated, the Nt triangle is assigned as charge-drain
% triangle
%   star: 10x(Nt-1) [ ... [T0 T1 T2 T3 ed1 ed2 ed3 v1 v2 v3].' ... ]
%   tri_star: it relates the triangle with the star basis function that
%   affect it 
%   tri_star: [ ... [ star_cent star_adj_1 i_star_ad_1 star_adj_2 i_star_adj_2 star_adj_3 i_star_adj_3 ].'  ... ]
%   i_* local index related to the central triangle of the star_adj_* basis function that affects
%   tri_star
%
% by Ed. Ubeda, nov. 2008

function [star, tri_star] = get_star_d3(obj);

Nt = length( obj.topol );

star = zeros(10,Nt-1);
tri_star = zeros(7,Nt);

for n=1:Nt-1,
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%% Star update %%%%
    %%%%%%%%%%%%%%%%%%%%%
    ed = obj.trian(1:3,n);
    i_nz_ed = find(ed~=0);
    
    ed = ed(i_nz_ed);
    s_ed = (1 + sign(ed))/2 + 1;
    ed = abs(ed);
    
    T = obj.edges( (ed-1)*4 + s_ed );
    v = obj.edges( (ed-1)*4 + 2 + s_ed );
    
    star(1,n) = n;
    star(1+i_nz_ed,n) = T;
    star(4+i_nz_ed,n) = ed;
    star(7+i_nz_ed,n) = v;
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    
    %% tri_star update
    %% central triangle
    tri_star(1,n) = n;
    
end;    %% for n=1:Nt-1,

%%% The basis functions must be assigned to the triangles taking into
%%% account their particular ordering 
for m=1:Nt,
    
    i_1 = find(star(2,:)==m);
    if length(i_1),
        for q=1:length(i_1),
            tri_star(1+2*(q-1)+1,m)=i_1(q);
            tri_star(1+2*(q-1)+2,m)=1;
        end;    %%% for q=length(i_1),
    end;    %% if length(i_1),
    
    i_2 = find(star(3,:)==m);
    if length(i_2),
        for q=1:length(i_2),
            tri_star(1+2*length(i_1)+2*(q-1)+1,m)=i_2(q);
            tri_star(1+2*length(i_1)+2*(q-1)+2,m)=2;
        end;    %%% for q=length(i_1),
    end;    %% if length(i_2),
    
    i_3 = find(star(4,:)==m);
    if length(i_3),
        for q=1:length(i_3),
            tri_star(1+2*(length(i_1)+length(i_2))+2*(q-1)+1,m)=i_3(q);
            tri_star(1+2*(length(i_1)+length(i_2))+2*(q-1)+2,m)=3;
        end;    %%% for q=length(i_1),    
    end;    %% if length(i_3),
    
end;  %%% for m=1:Nt,
