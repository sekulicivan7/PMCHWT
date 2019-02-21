% Script to compute the quantities related with the distribution of triangular
% prisms attached inside the surface triangulation
%
%  obj_inprismV = obj_inpV2
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
%    .un_nrm_v : [ 3 x Nv]
%      Surface-average OUTWARD unit vectors normal to the vertices in the surface
%
% by Ed. Ubeda, October 2013
%
%   .un_nrm_v is maintained
%   MOREOVER new fields are defined
%   .un_nrm_tri_v: [ 3 x [3xNt] ]
%    The OUTWARD normal vector depends on the vertex and the triangle
%    it arises from the average of the normal vectors to the adjacent triangles
%    to the edges defining that vertex
%   . cent_Spr: [3 x Nt]
%    Centroid of the bootom triangle
%
%  If (check_border_cross), the normal vectors are controlled so that they do not
%  cross the surface. For that they are compared against .un_nrm_v (which
%  establish the limit). If they do not overpass them, they are kept; if they
%  do, they are corrected.
%
% Modified by Ed. Ubeda, february 2014

function obj_inpV2 = get_obj_inprismV_Ivan2( obj , obj_V, loop_or , fact_vl )

check_border_cross = 1;

obj_inpV2 = struct('vertex',[],'topol_Qpr',[],'topol_tet_1',[],'topol_tet_2',[],'topol_tet_3',[],'topol_Qpr_T1',[],'topol_Qpr_T2',[],'topol_Spr',[],'ind_Sc',[],'un_int_Qpr_T1',[],'un_int_Qpr_T2',[],'un_int_Spr',[],'ds_Qpr_T1',[],'ds_Qpr_T2',[],'ds_Spr',[],'ln_Qpr',[],'fact',[],'Hpr',[],'un_nrm_v',[],'un_nrm_tri_v',[],'invrt_to_ed',[],'cent_Spr',[]);

obj_inpV2.fact = fact_vl;

Nt = length( obj.ds);
Ne = length( obj.edges );
Nv = length( obj.vertex );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (I) computation of the normal vector to each vertex by            %%
%%% Setting the surface-weighted average normal unit vector un_nrm_v  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_inpV2.un_nrm_v = zeros(3,Nv);
for m=1:Nv,
    %%% For each triangle sharing the m-th vertex
    %%% Compute the normal average triangle
    un_tmp = zeros(3,1);
    for q=1:loop_or.nv(m),
        Sq = loop_or.tri(m,q);
        %%%% un_tmp = un_tmp + obj.ds(Sq) * obj.un(:,Sq);
        %%%% REPLACED BY THE FOLLOWING DEFINITION !!!
        un_tmp = un_tmp + obj.un(:,Sq);
    end;  %% for q=1:loop_or.nv(m),
    
    obj_inpV2.un_nrm_v(:,m) = unitary(un_tmp);
    
end; %% for m=1:Nv,

%%% SETUP %%%
obj_inpV2.vertex = zeros(3,3*Nt);
obj_inpV2.Hpr = zeros(1,Nt);

obj_inpV2.topol_Qpr = zeros(4,3*Nt);
obj_inpV2.topol_Qpr_T1 = zeros(3,3*Nt);
obj_inpV2.topol_Qpr_T2 = zeros(3,3*Nt);
obj_inpV2.topol_Spr = zeros(3,Nt);

obj_inpV2.ind_Sc = zeros(2,Ne);
obj_inpV2.ln_Qpr = zeros(1,3*Nt);
obj_inpV2.invrt_to_ed = zeros(1,3*Nt);

obj_inpV2.topol_tet_1 = zeros(4,Nt);
obj_inpV2.topol_tet_2 = zeros(4,Nt);
obj_inpV2.topol_tet_3 = zeros(4,Nt);

%%%  new forget_obj_inprismV_2
%%% >>>
%%% >>>
obj_inpV2.un_nrm_tri_v = zeros( 3 , 3*Nt);
obj_inpV2.cent_Spr = zeros( 3 , Nt);
%%% <<<
%%% <<<

for p=1:Nt,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (II)  Computation of the set of vertices related to each     %%
    %% triangle and vertex of the surface meshing                   %%
    %% Computation of the average side length of the triangle       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ln = sum( obj.ln(abs(obj.trian(1:3,p))).' )./3;
    
    %%% computation of the height of the pentahedron in terms of Ln and
    %%% fact_vl
    Hn = fact_vl * Ln;
 
    
    %% update obj_inpV2.Hpr
    obj_inpV2.Hpr(p) = Hn;
    
    %%% OUTWARD normal vector to the triangle
    un_p = +obj.un(:,p);
    
    %%% Three vertices defining the pth surface triangle
    v_S = obj.topol(:,p);
    
    %% First, second and third vertex update on the pth triangle
    for n=1:3,
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% (III) Definition of the normal vector INWARDS along the intersecting segment         %%%
        %%%       between the two oposed triangles to the edges that share the particular vertex %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% new for get_obj_inprismV_2
        %%%% >>>
        %%%% >>>
        i_v_ed_S_mn_1 = mod( (n - 1) - 1, 3 ) + 1;
        i_v_ed_S_pl_1 = mod( (n - 1) + 1, 3 ) + 1;
        
        vert_or    = obj.vertex( : , obj.topol(n,p) );
        vert_pl_1 = obj.vertex( : , obj.topol( i_v_ed_S_pl_1 , p ) );
        vert_mn_1 = obj.vertex( : , obj.topol( i_v_ed_S_mn_1 , p ) );
        
        l_vec_pl_1 = unitary( vert_pl_1 - vert_or );
        l_vec_mn_1 = unitary( vert_mn_1 - vert_or );
        
        %% Minus-1 Adjacent triangle (note that surface minus1 corresponds with lvec plus)
        ed_mn_1 = obj.trian( i_v_ed_S_mn_1 , p );
        if ( sign(ed_mn_1)== +1 ),
            Tad_mn_1 = obj.edges( 2 , abs(ed_mn_1) );
        elseif ( sign(ed_mn_1)== -1 ),
            Tad_mn_1 = obj.edges( 1 , abs(ed_mn_1) );
        end;  %%% if ( sign(ed_mn_1) == +1 ),
        un_mn_1 = cross(  unitary( obj.un( : , Tad_mn_1 ) + obj.un( : , p ) ) , l_vec_pl_1 );
        
        %% plus 1 Adjacent triangle (note that surface plus1 corresponds with lvec minus1)
        ed_pl_1 = obj.trian( i_v_ed_S_pl_1 , p );
        if ( sign(ed_pl_1)== +1 ),
            Tad_pl_1 = obj.edges( 2 , abs(ed_pl_1) );
        elseif ( sign(ed_pl_1)== -1 ),
            Tad_pl_1 = obj.edges( 1 , abs(ed_pl_1) );
        end;  %%% if ( sign(ed_mn_1) == +1 ),
        un_pl_1 = cross( l_vec_mn_1 , unitary( obj.un( : , Tad_pl_1 ) + obj.un( : , p ) ) );
        
        vec_un_2 = unitary( un_pl_1 + un_mn_1 );
        vec_un_1 = unitary( un_pl_1 - un_mn_1 ) ;
        
        vec_nrm = cross( vec_un_2 , vec_un_1 );
        
        obj_inpV2.un_nrm_tri_v(:, 3*(p-1) + n ) = vec_nrm ;
        %%%% <<<
        %%%% <<<
        
        %%%% cos_ang_p_n = sum( un_p .* obj_inpV2.un_nrm_v(:,v_S(n) ));
        %%%% modified for get_obj_inprismV_2
        %%%% >>>
        cos_ang_p_n = sum( un_p .* vec_nrm );
        %%%% <<<
        H_p_n = Hn./cos_ang_p_n;
        
        %%%% vert_p_n = obj.vertex(:,v_S(n)) - H_p_n * obj_inpV2.un_nrm_v(:,v_S(n) );
        %%%% modified for get_obj_inprismV_2
        %%%% >>>
        vert_p_n = obj.vertex(:,v_S(n)) + H_p_n * vec_nrm;
        %%%% <<<
        
        %%% update n-th vertex of the p-th triangle
        obj_inpV2.vertex(:,3*(p-1) + n) = vert_p_n;
        
        %%% Find and Update the distance between the top vertex with the bottom vertex
        obj_inpV2.ln_Qpr(3*(p-1) + n) = norm( vert_p_n - obj.vertex(:,v_S(n)) );
        
        %%% Update the centroid of the bottom triangle
        obj_inpV2.cent_Spr(:,p) = obj_inpV2.cent_Spr(:,p) + vert_p_n/3;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% (IV)  Computation of the indices of the (quasi-)NORMAL           %%
        %%     Quadrilaterals (_Qpr) or of the First or Second Triangles    %%
        %%     (Qpr_T1, Qpr_T2) associated to each quadrilateral            %%
        %%     sharing a particular edge arising from the discretization.   %%
        %%     The row 1 refers to the Qpr or (Qpr_T1, Qpr_T2) linked to    %%
        %%     the T_plus of the edge.                                      %%
        %%     The row 2 refers to the Qpr or (Qpr_T1, Qpr_T2) linked to    %%
        %%     the T_min of the edge.                                       %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ed_n = obj.trian(n,p);
        s_ed_n = sign(ed_n);
        ed_n = abs(ed_n);
        
        if (s_ed_n==+1),
            i_n=1;
        elseif (s_ed_n==-1),
            i_n=2;
        end; %% if (s_ed_n==+1),  elseif (s_ed_n==-1),
        
        obj_inpV2.ind_Sc(i_n,ed_n) = 3*(p-1) + n;
        obj_inpV2.invrt_to_ed(3*(p-1) + n) = ed_n;
        
    end; %%% for n=1:3,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (V)     update of the vertex information in obj_inpV2.topol_Qpr,  %%
    %% obj_inpV2.topol_Qpr_T1, obj_inpV2.topol_Qpr_T2, obj_inpV2.topol_Spr %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n0 = (3*(p-1)); %% constant to take care of the vertices involved with the pth triangle
    
    %%% (a) Prism recangular faces, (pseudo-)NORMAL to the surface
    obj_inpV2.topol_Qpr(:,3*(p-1) + 1) =  [  v_S(2) -(n0+2)  -(n0+3)  v_S(3) ].';
    obj_inpV2.topol_Qpr(:,3*(p-1) + 2) =  [  v_S(3) -(n0+3)  -(n0+1)  v_S(1) ].';
    obj_inpV2.topol_Qpr(:,3*(p-1) + 3) =  [  v_S(1) -(n0+1)  -(n0+2)  v_S(2) ].';
    
    %%% (b) Triangle T1 of the prism rectangular faces, (pseudo-)NORMAL to
    %%% the surface
    obj_inpV2.topol_Qpr_T1(:,3*(p-1) + 1) = [ v_S(2) -(n0+2)  v_S(3) ].';
    obj_inpV2.topol_Qpr_T1(:,3*(p-1) + 2) = [ v_S(3) -(n0+3)  v_S(1) ].';
    obj_inpV2.topol_Qpr_T1(:,3*(p-1) + 3) = [ v_S(1) -(n0+1)  v_S(2) ].';
    
    %%% (c) Triangle T2 of the prism rectangular faces, (pseudo-)NORMAL to
    %%% the surface
    obj_inpV2.topol_Qpr_T2(:,3*(p-1) + 1) = [ -(n0+2)  -(n0+3)  v_S(3) ].';
    obj_inpV2.topol_Qpr_T2(:,3*(p-1) + 2) = [ -(n0+3)  -(n0+1)  v_S(1) ].';
    obj_inpV2.topol_Qpr_T2(:,3*(p-1) + 3) = [ -(n0+1)  -(n0+2)  v_S(2) ].';
    
    %%% (d) Prism basis
    obj_inpV2.topol_Spr(:,p) = [ -(n0+1) -(n0+2) -(n0+3) ].';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% (VI) Define the THREE tetrahedra that compose the prism (4 vertices involved for each) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj_inpV2.topol_tet_1(:,p) = [ v_S(1)  v_S(2)  v_S(3) -(n0+1) ].';
    obj_inpV2.topol_tet_2(:,p) = [ v_S(2)  v_S(3) -(n0+1) -(n0+2) ].';
    obj_inpV2.topol_tet_3(:,p) = [ v_S(3) -(n0+1) -(n0+2) -(n0+3) ].';
    
end; %%% for p=1:Nt,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (VII) Check whether the normal vectors in obj_inpV2.un_nrm_tri_v are  %%%%
%%% inside the boundary premises of the object;                          %%%%
%%% If they are not, the wedge the three bottom vertices of the wedge    %%%%
%%% are recomputed so that they fit in the object premises               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% new for get_obj_inprismV_2
%%%% >>>
%%%% >>>
if check_border_cross,
    
    for p=1:Nt,
        
        Hp = obj_inpV2.Hpr(p);
        v_p_S = obj.topol(:,p);
        cent_p = obj.cent(:,p);
        un_p = obj.un(:,p);
        
        for n=1:3,
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Each vertex check %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%
            cent_p_bottom = obj_inpV2.cent_Spr(:,p);
            un_p_n = obj_inpV2.un_nrm_tri_v( : , 3*(p-1) + n );  %%% (I)
            
            %%% Establish the vector_coordinates of the mid-plane
            %%% tt_vec = unitary( obj.vertex(:,v_p_S(n)) - cent_p );
            %%% nn_vec = cross( tt_vec , un_p );
            tt_1_vec = unitary( obj_inpV2.vertex( : , 3*(p-1) + n ) -  cent_p );
            tt_2_vec = unitary( obj.vertex(:,v_p_S(n)) -  cent_p_bottom );
            
            nn_vec = unitary( cross(tt_1_vec, tt_2_vec ) );
            
            %% projection on the plane of obj_inpV2.un_nrm_v(:,v_p_S(n))
            d_un_nrm_prj = sum( obj_inpV2.un_nrm_v(:,v_p_S(n)) .* nn_vec );
            un_nrm_v_prj = unitary( obj_inpV2.un_nrm_v(:,v_p_S(n)) - d_un_nrm_prj * nn_vec );
            
            %% projection on the plane of un_p
            d_p_prj = sum( un_p .* nn_vec );
            un_p_prj = unitary( un_p - d_p_prj * nn_vec );
            
            %%% compute the cosinus of the angles
            cos_ang_un_v = sum( un_p .* un_nrm_v_prj );
            cos_ang_ref =  sum( un_p .* un_p_n );
            
            if ( cos_ang_ref > cos_ang_un_v ),
                %%% WEDGE breaking into the boundary
                %%% The normal vectors in obj_inpV2.un_nrm_tri_v( : , 3*(p-1) + [1:3] )
                %%% need do be corrected
                
                %%% [a] Correction of the vertex directly affected
                if (cos_ang_un_v > +eps),
                    %% the un_nrm_v vector does not lie in the triangle surface or crosses the surface with less than 90 degrees angle
                    obj_inpV2.un_nrm_tri_v( : , 3*(p-1) + n ) = un_nrm_v_prj; %%% (II)
                elseif (cos_ang_un_v <= +eps),
                    %% Extreme case; the un_nrm_v vector lies in the triangle surface or it even goes down
                    %% The bottom triangle cannot be defined in this case under these premises
                    %% Therefore, the un_nrm_v_prj is redefined as the average of un_p_n and un_nrm_v_prj
                    %% ( even though it may cross the boundary a bit)
                    disp( 'Attention!!!: cos_ang_un_v smaller than zero!!!!' );
                    %%% I HAVE REALIZED THAT
                    %%% AFTHER REPLACING THE DEFINITION un_tmp = un_tmp + obj.ds(Sq) * obj.un(:,Sq);
                    %%% BY THE DEFINITION un_tmp = un_tmp + obj.un(:,Sq) AT
                    %%% THE NORMAL VECTORS TO THE VERTICES THIS SCENARIO
                    %%% (cos_ang_un_v <=0 ) IS NOT REACHED
                    %%%
                    un_nrm_v_prj = unitary( un_nrm_v_prj + un_p_n );
                    cos_ang_un_v = sum( un_p .* un_nrm_v_prj );
                end; %% if (cos_ang_un_v > +eps), elseif (cos_ang_un_v <= +eps),
                H_p_n = Hp./cos_ang_un_v;
                
                vert_p_n_0 = obj_inpV2.vertex( : , 3*(p-1) + n );
                
                vert_p_n = obj.vertex(:,v_p_S(n)) + H_p_n * un_nrm_v_prj;
                obj_inpV2.vertex( : , 3*(p-1) + n ) = vert_p_n;
                obj_inpV2.ln_Qpr( 3*(p-1) + n ) = norm( vert_p_n - obj.vertex(:,v_p_S(n)) );
                
                %%% [b] Correction of the other two vertices
                %%% Compute the centroid in the bottom triangle and the
                %%% reference values in the just corrected vertex
                l_vec_ref = vert_p_n - cent_p_bottom;
                l_nrm_ref = norm( l_vec_ref );
                
                l_vec_ref_0 = vert_p_n_0 - cent_p_bottom;
                l_nrm_ref_0 = norm( l_vec_ref_0 );
                
                %% scale factor
                scl_n  = l_nrm_ref./l_nrm_ref_0;
                
                %%% [n+1] vertex - correction
                np1 = mod( (n-1) + 1 , 3 ) + 1;
                vert_p_np1_0 = obj_inpV2.vertex( : , 3*(p-1) + np1 );
                l_vec_or_np1 = vert_p_np1_0 - cent_p_bottom;
                vert_p_np1 = cent_p_bottom  +  scl_n * l_vec_or_np1;
                obj_inpV2.vertex( : , 3*(p-1) + np1 ) = vert_p_np1;
                vec_tmp_np1 = obj.vertex(:,v_p_S(np1)) - vert_p_np1;
                obj_inpV2.un_nrm_tri_v( : , 3*(p-1) + np1 ) = unitary( vec_tmp_np1 );  %%% (III)
                obj_inpV2.ln_Qpr(3*(p-1) + np1) = norm( vec_tmp_np1 );
                
                %%% [n-1] vertex - correction
                nm1 = mod( (n-1) - 1 , 3 ) + 1;
                vert_p_nm1_0 = obj_inpV2.vertex( : , 3*(p-1) + nm1 );
                l_vec_or_nm1 = vert_p_nm1_0 - cent_p_bottom;
                vert_p_nm1 = cent_p_bottom  +  scl_n * l_vec_or_nm1;
                obj_inpV2.vertex( : , 3*(p-1) + nm1 ) = vert_p_nm1;
                vec_tmp_nm1 = obj.vertex(:,v_p_S(nm1)) - vert_p_nm1;
                obj_inpV2.un_nrm_tri_v( : , 3*(p-1) + nm1 ) = unitary( vec_tmp_nm1 );   %%% (IV)
                obj_inpV2.ln_Qpr(3*(p-1) + nm1) = norm( vec_tmp_nm1 );
                
                %%% update the centroid at the bottom triangle
                obj_inpV2.cent_Spr(:,p) =  ( vert_p_n + vert_p_np1 + vert_p_nm1 )./ 3;
                
            end; %% if ( cos_ang_ref < cos_ang_un_v ),
            
        end; %%% for n=1:3,
        
    end; %%% for p=1:Nt,
    
end; %%% if check_border_cross,
%%%% <<<
%%%% <<<

%% Setups
obj_inpV2.un_int_Qpr_T1 = zeros(3,3*Nt);
obj_inpV2.un_int_Qpr_T2 = zeros(3,3*Nt);
obj_inpV2.un_int_Spr = zeros(3,Nt);
obj_inpV2.ds_Qpr_T1 = zeros(1,3*Nt);
obj_inpV2.ds_Qpr_T2 = zeros(1,3*Nt);
obj_inpV2.ds_Spr = zeros(1,Nt);

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
        v1 = obj.vertex( : , obj_inpV2.topol_Qpr_T1( 1 , 3*(p-1) + n ) );
        v2 = obj_inpV2.vertex( :  , abs( obj_inpV2.topol_Qpr_T1( 2 , 3*(p-1) + n ) ) );
        v3 = obj.vertex( : , obj_inpV2.topol_Qpr_T1( 3, 3*(p-1) + n ) );
        c    = cross(v3-v1,v2-v1);
        obj_inpV2.un_int_Qpr_T1( : , 3*(p-1) + n) = -unitary(c);% MINUS ADDED, NORMALS OF THE SIDE QUADRILATERALS NEED TO BE IN THE SAME DIRECTION AS WHEN THEY ARE IN
                                                                    % THE OBJECT
        obj_inpV2.ds_Qpr_T1(3*(p-1) + n)   = sqrt(sum(c.^2))/2;
    end; %%% for n=1:3,
    
    %%% (b) Triangle T2 of the (pseudo-)NORMAL quadrilaterals
    for n=1:3,
        v1 = obj_inpV2.vertex( : , abs(obj_inpV2.topol_Qpr_T2( 1 , 3*(p-1) + n )) );
        v2 = obj_inpV2.vertex( :  , abs(obj_inpV2.topol_Qpr_T2(2,3*(p-1) + n)) );
        v3 = obj.vertex( : , obj_inpV2.topol_Qpr_T2( 3, 3*(p-1) + n ) );
        c    = cross(v3-v1,v2-v1);
        obj_inpV2.un_int_Qpr_T2( : , 3*(p-1) + n) = -unitary(c);% MINUS ADDED, NORMALS OF THE SIDE QUADRILATERALS NEED TO BE IN THE SAME DIRECTION AS WHEN THEY ARE IN
                                                                    % THE OBJECT
        obj_inpV2.ds_Qpr_T2(3*(p-1) + n)   = sqrt(sum(c.^2))/2;
    end; %%% for n=1:3,

    %%% (c) Opposed Triangle (parallel) to the surface triangle
    v1 = obj_inpV2.vertex( : , abs( obj_inpV2.topol_Spr( 1 , p) ) );
    v2 = obj_inpV2.vertex( :  , abs( obj_inpV2.topol_Spr( 2 ,p) ) );
    v3 = obj_inpV2.vertex( : , abs( obj_inpV2.topol_Spr( 3, p) ) );
    c    = cross( v3-v1 , v2-v1 );
    obj_inpV2.un_int_Spr( : , p ) = unitary( c );
    obj_inpV2.ds_Spr(p)   = sqrt(sum(c.^2))/2;
    
end

