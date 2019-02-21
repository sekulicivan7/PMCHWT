%
% Function that provides the unit tangential vectors and the origin vertices along the edges
%
% l_vec = 3xNe
%   - unitary vector along the edge (along the rotation sense)
%  vert_or = 3xNe
%   - original vertex of the edges
%
% by Ed. Ubeda, october 2012

function [ l_vec, vert_or, vert_mid, un_c_plus, un_c_min] = get_l_vec(obj)

Ne = length(obj.ln);

l_vec = [];
vert_or = [];
vert_mid = [];

un_c_plus = [];
un_c_min = [];

for m=1:Ne,
    
    Tp = obj.edges(1,m);
    
    v_op_p = obj.edges(3,m);
    v_all_p = obj.topol(:,Tp);
    i_op = find(v_all_p==v_op_p);
    %%% rotation sense (first: +1 ; second: +2
    i_common = [ rem(i_op+1-1,3)+1 rem(i_op+2-1,3)+1 ];
    v_common = v_all_p(i_common);
    
    l_tmp_vec = unitary( obj.vertex(:,v_common(2)) - obj.vertex(:,v_common(1)) );
    l_vec = [ l_vec  l_tmp_vec ];
    
    vert_or = [ vert_or obj.vertex(:,v_common(1)) ];
    vert_mid = [ vert_mid   0.5* ( obj.vertex(:,v_common(2)) + obj.vertex(:,v_common(1)) ) ];
    
    un_c_tmp_vec_plus = cross( obj.un(:,Tp) , l_tmp_vec );
    un_c_plus = [ un_c_plus un_c_tmp_vec_plus ];

    Tm = obj.edges(2,m);    
    un_c_tmp_vec_min = cross( obj.un(:,Tm) , -l_tmp_vec );    
    un_c_min = [ un_c_min un_c_tmp_vec_min ];    
    
end;  %%% for m=1:Ne,

