% Script to compute the quantities related with the distribution of triangular
% prisms attached inside the surface triangulation
%
%  obj_inprismV = obj_nonconf_inpV
%
%   . vertex: [ 3 x (3xNt) ]
%      List of vertices inside the volume (aggregation of three tetrahedra)
%   . topol_tet_1 : [ 4 x Nt ]
%      List of 4 vertices for the FIRST tetrahedron
%     [ ... [ +v1; +v2; +v3; -v4 ] ... ]
%   . topol_tet_2 : [ 4 x Nt ]
%      List of 4 vertices for the SECOND tetrahedron
%     [ ... [ +v2; +v3; -v4; -v5 ] ... ]
%   . topol_tet_3 : [ 4 x Nt ]
%      List of 4 vertices for the THIRD tetrahedron
%     [ ... [ +v3; -v4; -v5; -v6 ] ... ]
%   . topol_Qpr: [4 x (3xNt)]
%      set of 4 vertices for each (pseudo-)NORMAL quadrilateral in the prism (3xNt)
%       [ ... [+v1;+v2;-v3;-v4] ... ]
%           +v1,+v2: indices of obj.vertex (range +(1:Nv))
%           -v3,-v4: indices of obj_inprismV.vertex  (range -(1:3*Nt))
%           [ The rotation convention is INWARDS of the prism ]
%   . topol_Qpr_T1: [3 x (3xNt)]
%      set of 3 vertices of the FIRST triangle for each (pseudo-)NORMAL quadrilateral in the prism (3xNt)
%       [ ... [+v1;+v2;-v3] ... ]
%           +v1,+v2: indices of obj.vertex  (range +(1:Nv))
%           -v3: indices of obj_inprismV.vertex (range -(1:3*Nt))
%           [ The rotation convention is INWARDS of the prism]
%   . topol_Qpr_T2: [3 x (3xNt)]
%      set of 3 vertices of the SECOND triangle for each (pseudo-)NORMAL quadrilateral in the prism (3xNt)
%       [ ... [-v3;-v4;+v1] ... ]
%           -v3,-v4: indices of obj_inprismV.vertex (range -(1:3*Nt)) 
%           +v1: indices of obj.vertex (range +(1:Nv))
%           [ The rotation convention is INWARDS of the prism]
%   . topol_Spr: [3 x Nt]
%      set of 3 vertices for each opposed basis prism triangle (PARALLEL to the surface)
%       (indices of vert_Qpr)
%       [ ... [-va;-vb;-vc] ... ] (range -(1:3*Nt)) 
%       [ The rotation convention is INWARDS of the prism]
%   . ind_Sc: [ 2 x Ne ]
%       [ ... [ iQ1;iQ2 ] ... ]
%      Two quadrilateral surfaces, at each side of the common edge
%       iQ1 and iQ2 are indices of the three quadriteral surfaces Qpr (or corresponding triangles Qpr_T1 - Qpr_T2)  
%       for each surface triangle: The range of iQ1-iQ2 is [1:3xNt]
%         - The row 1 refers to the Qpr or (Qpr_T1, Qpr_T2) associated with the T_plus of the edge in obj.edges
%         - The row 2 refers to the Qpr or (Qpr_T1, Qpr_T2) associated with the T_minus of the edge in obj.edges
%   . invrt_to_ed: [1 x (3xNt)]
%       It relates the interior-prism index of the triangle with the oposed edge
%   . un_int_Qpr_T1: [ 3 x (3xNt) ]
%      Outer normal vector to the FIRST triangle of the (Pseudo-)NORMAL quadrilateral around the prism
%   . un_int_Qpr_T2: [ 3 x (3xNt) ]
%      Outer normal vector to the SECOND triangle of the (Pseudo-)NORMAL quadrilateral around the prism
%   . un_int_Spr: [ 3 x Nt ]
%      Outer normal vector to the triangle opposed to the surface triangle
%   . ds_Qpr_T1: [ 1 x (3xNt) ]
%      Surface of each FIRST triangle of each (pseudo-)NORMAL quadilateral in each prism
%   . ds_Qpr_T2: [ 1 x (3xNt) ]
%      Surface of each SECOND triangle of each (pseudo-)NORMAL quadilateral in each prism
%   . ds_Spr: [ 1 x Nt ]
%      Surface of each oposed basis prism triangle
%   . ln_Qpr: [ 1 x (3xNt) ] 
%      length of each segment connecting each vertex of the top surface triangle 
%      with the the bottom triangle
%   . fact: number
%      factor of the average length of the 3 edges of each surface triangle
%      that produces the height of the prisms
%   . Hpr : [ 1 x Nt]
%      Height of each prism
%
%   Modified by Ed. Ubeda, February 2014
%   +  The inward vector is perpendicular to the surface 
%      (the wedges become rectangular wedges)%
%   + The information u_nc_nrm is ignored
%
% by Ed. Ubeda, October 2013
%
% [ASSUMED right triangular prisms]
% 

function obj_nc_inpV = get_obj_nc_inprismV_Ivan( obj , fact_vol );

%%% obj_nc_inpV = struct('vertex',[],'topol_Qpr',[],'topol_tet_1',[],'topol_tet_2',[],'topol_tet_3',[],'topol_Qpr_T1',[],'topol_Qpr_T2',[],'topol_Spr',[],'ind_Sc',[],'un_int_Qpr_T1',[],'un_int_Qpr_T2',[],'un_int_Spr',[],'ds_Qpr_T1',[],'ds_Qpr_T2',[],'ds_Spr',[],'ln_Qpr',[],'fact',[],'Hpr',[],'un_nrm_v',[],'invrt_to_ed',[]);
obj_nc_inpV = struct('vertex',[],'topol_Qpr',[],'topol_tet_1',[],'topol_tet_2',[],'topol_tet_3',[],'topol_Qpr_T1',[],'topol_Qpr_T2',[],'topol_Spr',[],'ind_Sc',[],'un_int_Qpr_T1',[],'un_int_Qpr_T2',[],'un_int_Spr',[],'ds_Qpr_T1',[],'ds_Qpr_T2',[],'ds_Spr',[],'ln_Qpr',[],'fact',[],'Hpr',[],'un_nrm_tri_v',[]);

obj_nc_inpV.fact = fact_vol;

Nt = length( obj.ds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (I) computation of the normal vector to each vertex by            %%
%%% Setting the surface-weighted average normal unit vector un_nrm_v  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % obj_nc_inpV.un_nrm_v = zeros(3,Nv);
% % 
% % for m=1:Nv,
% %     
% %     %%% For each triangle sharing the m-th vertex
% %     %%% Compute the normal average triangle
% %     un_tmp = zeros(3,1);
% %     for q=1:loop_or.nv(m),
% %         Sq = loop_or.tri(m,q);
% %         un_tmp = un_tmp + obj.ds(Sq) * obj.un(:,Sq);
% %     end;  %% for q=1:loop_or.nv(m),
% %     
% %     obj_nc_inpV.un_nrm_v(:,m) = unitary(un_tmp);
% %     
% % end; %% for m=1:Nv,


%%% SETUP %%%
obj_nc_inpV.vertex = zeros(3,3*Nt);
obj_nc_inpV.Hpr = zeros(1,Nt);

obj_nc_inpV.topol_Qpr = zeros(4,3*Nt);
obj_nc_inpV.topol_Qpr_T1 = zeros(3,3*Nt);
obj_nc_inpV.topol_Qpr_T2 = zeros(3,3*Nt);
obj_nc_inpV.topol_Spr = zeros(3,Nt);

% % obj_nc_inpV.ind_Sc = zeros(2,Ne);
obj_nc_inpV.ln_Qpr = zeros(1,3*Nt);
obj_nc_inpV.invrt_to_ed = zeros(1,3*Nt);

obj_nc_inpV.topol_tet_1 = zeros(4,Nt);
obj_nc_inpV.topol_tet_2 = zeros(4,Nt);
obj_nc_inpV.topol_tet_3 = zeros(4,Nt);

L1 = sqrt( sum( ( obj.vertex(:,obj.topol(3,:)) - obj.vertex(:,obj.topol(2,:)) ).^2 ) );
L2 = sqrt( sum( ( obj.vertex(:,obj.topol(3,:)) - obj.vertex(:,obj.topol(1,:)) ).^2 ) );
L3 = sqrt( sum( ( obj.vertex(:,obj.topol(2,:)) - obj.vertex(:,obj.topol(1,:)) ).^2 ) );

obj_nc_inpV.cent_Spr = zeros( 3 , Nt);
obj_nc_inpV.un_nrm_tri_v = zeros( 3 , 3*Nt);

for p=1:Nt,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (II)  Computation of the set of vertices related to each     %%
    %% triangle and vertex of the surface meshing                   %%
    %% Computation of the average side length of the triangle       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     Ln = sum( obj.ln(abs(obj.trian(1:3,p))).' )./3;
  Ln = ( L1(p) + L2(p) + L3(p) )./3;
    
    %%% computation of the height of the pentahedron in terms of Ln and
    %%% fact_vol
    Hp = fact_vol * Ln;
    
    %% update obj_nc_inpV.Hpr
    obj_nc_inpV.Hpr(p) = Hp;
    
    %%% OUTWARD normal vector to the triangle
    un_p = +obj.un(:,p);
    
    %%% Three vertices defining the pth surface triangle
    v_S = obj.topol(:,p);
    
    %% First, second and third vertex update on the pth triangle
    for n=1:3,
        
        obj_nc_inpV.un_nrm_tri_v(:, 3*(p-1) + n ) = un_p ;  
        
        vert_p_n = obj.vertex(:,v_S(n)) - Hp * un_p;
                
        %%% update n-th vertex of the p-th triangle
        obj_nc_inpV.vertex(:,3*(p-1) + n) = vert_p_n;
        
        %%% Find and Update the distance between the top vertex with the bottom vertex
        obj_nc_inpV.ln_Qpr(3*(p-1) + n) = norm( vert_p_n - obj.vertex(:,v_S(n)) );
        
        %%% Update the centroid of the bottom triangle
        obj_nc_inpV.cent_Spr(:,p) = obj_nc_inpV.cent_Spr(:,p) + vert_p_n/3;
        
        %%% OUT FOR get_obj_nc_inprismV
        %%% >>>
        %%% >>>        
        
        % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % %         %% (IV)  Computation of the indices of the (quasi-)NORMAL           %%
        % % %         %%     Quadrilaterals (_Qpr) or of the First or Second Triangles    %%
        % % %         %%     (Qpr_T1, Qpr_T2) associated to each quadrilateral            %%
        % % %         %%     sharing a particular edge arising from the discretization.   %%
        % % %         %%     The row 1 refers to the Qpr or (Qpr_T1, Qpr_T2) linked to    %%
        % % %         %%     the T_plus of the edge.                                      %%
        % % %         %%     The row 2 refers to the Qpr or (Qpr_T1, Qpr_T2) linked to    %%
        % % %         %%     the T_min of the edge.                                       %%
        % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % %         ed_n = obj.trian(n,p);
        % % %         s_ed_n = sign(ed_n);
        % % %         ed_n = abs(ed_n);
        % % %         
        % % %         if (s_ed_n==+1),
        % % %             i_n=1;
        % % %         elseif (s_ed_n==-1),
        % % %             i_n=2;
        % % %         end; %% if (s_ed_n==+1),  elseif (s_ed_n==-1),
        % % %         
        % % %         obj_nc_inpV.ind_Sc(i_n,ed_n) = 3*(p-1) + n;
        % % %         obj_nc_inpV.invrt_to_ed(3*(p-1) + n) = ed_n;
        
        %%% <<<
        %%% <<<
        
    end; %%% for n=1:3,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (V)     update of the vertex information in obj_nc_inpV.topol_Qpr,  %%
    %% obj_nc_inpV.topol_Qpr_T1, obj_nc_inpV.topol_Qpr_T2, obj_nc_inpV.topol_Spr %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    n0 = (3*(p-1)); %% constant to take care of the vertices involved with the pth triangle

    %%% (a) Prism recangular faces, (pseudo-)NORMAL to the surface
    obj_nc_inpV.topol_Qpr(:,3*(p-1) + 1) =  [  v_S(2) -(n0+2)  -(n0+3)  v_S(3) ].';
    obj_nc_inpV.topol_Qpr(:,3*(p-1) + 2) =  [  v_S(3) -(n0+3)  -(n0+1)  v_S(1) ].';
    obj_nc_inpV.topol_Qpr(:,3*(p-1) + 3) =  [  v_S(1) -(n0+1)  -(n0+2)  v_S(2) ].';
    
    %%% (b) Triangle T1 of the prism rectangular faces, (pseudo-)NORMAL to
    %%% the surface
    obj_nc_inpV.topol_Qpr_T1(:,3*(p-1) + 1) = [ v_S(2) -(n0+2)  v_S(3) ].';
    obj_nc_inpV.topol_Qpr_T1(:,3*(p-1) + 2) = [ v_S(3) -(n0+3)  v_S(1) ].';
    obj_nc_inpV.topol_Qpr_T1(:,3*(p-1) + 3) = [ v_S(1) -(n0+1)  v_S(2) ].';
    
    %%% (c) Triangle T2 of the prism rectangular faces, (pseudo-)NORMAL to
    %%% the surface    
    obj_nc_inpV.topol_Qpr_T2(:,3*(p-1) + 1) = [ -(n0+2)  -(n0+3)  v_S(3) ].';
    obj_nc_inpV.topol_Qpr_T2(:,3*(p-1) + 2) = [ -(n0+3)  -(n0+1)  v_S(1) ].';
    obj_nc_inpV.topol_Qpr_T2(:,3*(p-1) + 3) = [ -(n0+1)  -(n0+2)  v_S(2) ].';
    
    %%% (d) Prism basis
    obj_nc_inpV.topol_Spr(:,p) = [ -(n0+1) -(n0+2) -(n0+3) ].'; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% (VI) Define the THREE tetrahedra that compose the prism (4 vertices involved for each) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj_nc_inpV.topol_tet_1(:,p) = [ v_S(1)  v_S(2)  v_S(3) -(n0+1) ].';
    obj_nc_inpV.topol_tet_2(:,p) = [ v_S(2)  v_S(3) -(n0+1) -(n0+2) ].';
    obj_nc_inpV.topol_tet_3(:,p) = [ v_S(3) -(n0+1) -(n0+2) -(n0+3) ].';
    
end; %%% for p=1:Nt,

%% Setups
obj_nc_inpV.un_int_Qpr_T1 = zeros(3,3*Nt);
obj_nc_inpV.un_int_Qpr_T2 = zeros(3,3*Nt);
obj_nc_inpV.un_int_Spr = zeros(3,Nt);
obj_nc_inpV.ds_Qpr_T1 = zeros(1,3*Nt);
obj_nc_inpV.ds_Qpr_T2 = zeros(1,3*Nt);
obj_nc_inpV.ds_Spr = zeros(1,Nt);

for p=1:Nt,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (VIII) Computation of the OUTWARD normal vectors and of the         %%
    %%      triangle areas of the triangles (opposed triangles and T1,T2   %%
    %%      of the (pseudo-)normal quadrilaterals                          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % <<< reminder from get_edge3_linear_2 >>>
    % c    = cross(v3-v1,v2-v1);
    % obj.un   = unitary(c);
    % obj.ds   = sqrt(sum(c.^2))/2;
    % <<< -------------------------------- >>>    
    %%% (a) Triangle T1 of the (pseudo-)NORMAL quadrilaterals
    for n=1:3,
        v1 = obj.vertex( : , obj_nc_inpV.topol_Qpr_T1( 1 , 3*(p-1) + n ) );
        v2 = obj_nc_inpV.vertex( :  , abs( obj_nc_inpV.topol_Qpr_T1( 2 , 3*(p-1) + n ) ) );
        v3 = obj.vertex( : , obj_nc_inpV.topol_Qpr_T1( 3, 3*(p-1) + n ) ); 
        c    = cross(v3-v1,v2-v1);
        obj_nc_inpV.un_int_Qpr_T1( : , 3*(p-1) + n) = unitary(c); 
        obj_nc_inpV.ds_Qpr_T1(3*(p-1) + n)   = sqrt(sum(c.^2))/2;
    end; %%% for n=1:3,
    
    %%% (b) Triangle T2 of the (pseudo-)NORMAL quadrilaterals
    for n=1:3,
        v1 = obj_nc_inpV.vertex( : , abs(obj_nc_inpV.topol_Qpr_T2( 1 , 3*(p-1) + n )) );
        v2 = obj_nc_inpV.vertex( :  , abs(obj_nc_inpV.topol_Qpr_T2(2,3*(p-1) + n)) );
        v3 = obj.vertex( : , obj_nc_inpV.topol_Qpr_T2( 3, 3*(p-1) + n ) ); 
        c    = cross(v3-v1,v2-v1);
        obj_nc_inpV.un_int_Qpr_T2( : , 3*(p-1) + n) = unitary(c);
        obj_nc_inpV.ds_Qpr_T2(3*(p-1) + n)   = sqrt(sum(c.^2))/2;
    end; %%% for n=1:3,
    
    %%% (c) Opposed Triangle (parallel) to the surface triangle
    v1 = obj_nc_inpV.vertex( : , abs( obj_nc_inpV.topol_Spr( 1 , p) ) );
    v2 = obj_nc_inpV.vertex( :  , abs( obj_nc_inpV.topol_Spr( 2 ,p) ) );
    v3 = obj_nc_inpV.vertex( : , abs( obj_nc_inpV.topol_Spr( 3, p) ) ); 
    c    = cross( v3-v1 , v2-v1 );
    obj_nc_inpV.un_int_Spr( : , p ) = unitary( c );
    obj_nc_inpV.ds_Spr(p)   = sqrt(sum(c.^2))/2;
    
end; %%% for p=1:Nt,