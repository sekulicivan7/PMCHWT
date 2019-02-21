%   Procedure to obtain the volumetric-tetrahedral parameters
%   associated to the interior volume
%   
%   + Input parameterss:
%       
%       vert_mid:  [3xNe] mid edge vertex
%       un_c_plus: [3xNe] outward contour normal vector from the positive side of the edge
%       un_c_min:  [3xNe] outward contour normal vector from the negative side of the edge
%       l_vec      [3xNe] tangential vector to the edge turning inwards
%       flag_dpth: [1xNe]  
%                  if 0, the depth of the interior vertex is established from the length of the common edge (by default: depth = ln/10)
%                  if not zero, it sets the value of the length of the depth of the interior edge with the common edge 
%
%   + Output parameters:
%
%       obj_V.
%           vert_ad  [ 3 x Ne ]
%              + Added vertex for each original edge in the original
%              surface mesh (INSIDE the object)
%           vert_out [3xNe]
%              + Added vertex for each original edge in the original
%              surface mesh (OUTSIDE the object; symmetric wrt the surface to vert_ad)
%           topol_ad [ 3 x (5xNe) ]
%              + Set of added 5 triangles for each basis function in the original surface mesh
%              [ the vertex rotation turns ouwards --> obj_v.un_ad ]
%           cent_ad  [ 3 x (5xNe) ]
%              + Centroids of the added triangles
%           ds_ad    [ 1 x (5xNe) ]
%              + Surface of the added triangles
%           un_ad    [ 3 x (5xNe) ]
%              + Outward vectors of the added triangles
%           monoSWG  [ 8 x (2xNe) ]
%              + Parameters defining the monoSWG basis functions
%              + if index > Ne;  corresponding to Tm in the original surface mesh
%              + if index <= Ne; corresponding to Tp in the original surface mesh
%              + [ (+)T1 (-)T2 (-)T3 (-)T4 (+-)v1 (+-)v2 (+-)v3 (+-)v4] .' * ones(1,2*Ne) 
%                   + --> triangle/vertex in the original mesh (obj)
%                   - --> triangle/vertex in the added mesh (obj_V)
%           Tmid    [ 1xNe ]
%               + Border triangle between the positive and the negative
%               tetrahedra
%           cent_V   [ 3 x [2xNe] ]
%               + Centroids of the 2xNe tetrahedra
%           dV       [ 1 x [2xNe] ]
%               + Volumes of the 2xNe tetrahedra
%           dV_nrm:  [ 1 x (2xNe) ]
%           dS_ad_nrm:  [ 1 x Ne ]
%               + Normalized surface of the central adjacent surface 
%   
%   by Ed. Ubeda, november 2012

function obj_V = get_obj_vol_Ivan( obj , un_c_plus , un_c_min , vert_mid , l_vec , fact_vl )

Ne = length( obj.ln );

obj_V = struct('vert_ad',[],'vert_out',[],'vert_S_pl',[],'vert_S_mn',[],'topol_ad',[],'monoSWG',[],'Tmid',[],'Tp_ad1',[],'Tp_ad2',[],'Tm_ad1',[],'Tm_ad2',[],'un_ed_int',[],'ds_ad',[],'depth_int',[],'cent_ad',[],'cent_V',[],'dV',[],'dV_nrm',[],'ds_ad_nrm',[],'vert_plIn',[],'vert_plOut',[],'vert_mnIn',[],'vert_mnOut',[],'vert_PL_pl',[],'vert_PL_out',[],'vert_PL_mn',[],'vert_PL_in',[],'vert_PL_plIn',[],'vert_PL_plOut',[],'vert_PL_mnIn',[],'vert_PL_mnOut',[],'vert_MN_plIn',[],'vert_MN_plOut',[],'vert_MN_mnIn',[],'vert_MN_mnOut',[],'vert_MN_pl',[],'vert_MN_out',[],'vert_MN_mn',[],'vert_MN_in',[]);

%%% Get the interior interior normal vectors
%%% obj_V.un_ed_int = unitary( cross( un_c_plus , l_vec ) + cross( l_vec , un_c_min ) );

% % % obj_V.un_ed_int = unitary( cross( un_c_plus , l_vec ) + cross( l_vec , un_c_min ) );
un_int_plus = cross( un_c_plus , l_vec ); %% inwards vector PL side of the edge
un_int_min = cross( l_vec , un_c_min );   %% inwards vector MN side of the edge
obj_V.un_ed_int = unitary( un_int_plus + un_int_min );

un_plOut = unitary( - un_c_plus +  obj_V.un_ed_int );
un_mnOut = unitary( - un_c_min +  obj_V.un_ed_int );
un_plIn = unitary( - un_c_plus -  obj_V.un_ed_int );
un_mnIn = unitary( - un_c_min -  obj_V.un_ed_int );

%%% Get the interior distribution of added vertices associated to the mid-points

depth_length = fact_vl*obj.ln;

%%%% TESTING VERTICES
%%%% >>>>
obj_V.vert_ad = vert_mid + ( ones(3,1) * depth_length ) .* obj_V.un_ed_int; 
obj_V.vert_out = vert_mid - ( ones(3,1) * depth_length ) .* obj_V.un_ed_int;
obj_V.vert_S_pl = vert_mid  - ( ones(3,1) * depth_length ) .* un_c_plus; 
obj_V.vert_S_mn = vert_mid - ( ones(3,1) * depth_length ) .* un_c_min;     
%
obj_V.vert_plOut = vert_mid + ( ones(3,1) * depth_length ) .* un_plOut; 
obj_V.vert_mnIn  = vert_mid + ( ones(3,1) * depth_length ) .* un_mnIn;     
obj_V.vert_mnOut = vert_mid + ( ones(3,1) * depth_length ) .* un_mnOut;        
obj_V.vert_plIn  = vert_mid + ( ones(3,1) * depth_length ) .* un_plIn;  
%
% DOUBLE THE 8 DIRECTIONS
%% [PL SIDE OF THE EDGE]
obj_V.vert_PL_pl = vert_mid + ( ones(3,1) * depth_length ) .* ( - un_c_plus );
obj_V.vert_PL_mn  = vert_mid + ( ones(3,1) * depth_length ) .* un_c_plus;
obj_V.vert_PL_in  = vert_mid + ( ones(3,1) * depth_length ) .* un_int_plus;
obj_V.vert_PL_out = vert_mid + ( ones(3,1) * depth_length ) .* ( - un_int_plus );
%
obj_V.vert_PL_plOut = vert_mid + ( ones(3,1) * depth_length ) .* unitary( - un_c_plus - un_int_plus  ); 
obj_V.vert_PL_mnIn  = vert_mid + ( ones(3,1) * depth_length ) .* unitary( un_c_plus +  un_int_plus );
obj_V.vert_PL_mnOut = vert_mid + ( ones(3,1) * depth_length ) .* unitary( un_c_plus -  un_int_plus );
obj_V.vert_PL_plIn  = vert_mid + ( ones(3,1) * depth_length ) .* unitary( - un_c_plus + un_int_plus  ); 
%% [MN SIDE OF THE EDGE]
obj_V.vert_MN_pl = vert_mid + ( ones(3,1) * depth_length ) .* un_c_min;
obj_V.vert_MN_mn  = vert_mid + ( ones(3,1) * depth_length ) .* ( - un_c_min );
obj_V.vert_MN_in  = vert_mid + ( ones(3,1) * depth_length ) .* un_int_min;
obj_V.vert_MN_out = vert_mid + ( ones(3,1) * depth_length ) .* ( - un_int_min );
%
obj_V.vert_MN_plOut = vert_mid + ( ones(3,1) * depth_length ) .* unitary( un_c_min - un_int_min );
obj_V.vert_MN_mnIn  = vert_mid + ( ones(3,1) * depth_length ) .* unitary( - un_c_min + un_int_min ); 
obj_V.vert_MN_mnOut = vert_mid + ( ones(3,1) * depth_length ) .* unitary( - un_c_min - un_int_min );
obj_V.vert_MN_plIn  = vert_mid + ( ones(3,1) * depth_length ) .* unitary(  un_c_min + un_int_min );  
%%%% <<<<

obj_V.depth_int = depth_length;

%% Get the set of added triangles / added normal vectors / added surfaces / added centroids
obj_V.topol_ad = zeros( 3 , 5*Ne );
obj_V.un_ad = zeros( 3 , 5*Ne );
obj_V.ds_ad = zeros( 1 , 5*Ne );
obj_V.ds_ad_nrm = zeros( 1 , Ne );
obj_V.Tmid = zeros( 1 , Ne );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    (I) ADDED-SURFACE PARAMETERS    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:Ne,
    
    %%% Positive Triangle - parameters
    Tp = obj.edges(1,m);
    vop_p = obj.edges(3,m);
    i_all_v_p = obj.topol(:,Tp);
    
    i_op_p = find( i_all_v_p == vop_p );
    
    vp_1 = obj.topol( rem( i_op_p + 1 - 1 , 3 ) + 1 , Tp );
    vp_2 = obj.topol( rem( i_op_p + 2 - 1 , 3 ) + 1 , Tp );
    
    %%% Negative Triangle - parameters
    Tm = obj.edges(2,m);
    vop_m = obj.edges(4,m);    
    i_all_v_m = obj.topol(:,Tm);
    
    i_op_m = find( i_all_v_m == vop_m );
    
    %%% swapped order in Vm wrt Tm to let the oposed triangles face-to-face int the corresponding vertices
    vm_1 = obj.topol( rem( i_op_m + 2 - 1 , 3 ) + 1 , Tm );
    vm_2 = obj.topol( rem( i_op_m + 1 - 1 , 3 ) + 1 , Tm );    
    
    %%% Update obj_V.topol
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% [First] Added-positive triangle  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    obj_V.topol_ad( : , (m-1)*5 + 1 ) = [ vop_p  vp_1  (-m) ];
    term_cross_1 = cross( ( obj.vertex(:,vp_1) - obj.vertex(:,vop_p) ) , ( obj_V.vert_ad(:,m) - obj.vertex(:,vop_p) ) );
    obj_V.un_ad( : , (m-1)*5 + 1 )  = unitary( term_cross_1 );   
    obj_V.ds_ad( (m-1)*5 + 1 ) = sqrt( sum(term_cross_1.^2) )./2;
    obj_V.cent_ad( : , (m-1)*5 + 1 ) = (1/3)*( obj.vertex(:,vop_p) + obj.vertex(:,vp_1) + obj_V.vert_ad(:,m) );        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% [Second] Added-positive triangle  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    obj_V.topol_ad( : , (m-1)*5 + 2 ) = [ vop_p  (-m)  vp_2 ];
    term_cross_2 = cross( ( obj_V.vert_ad(:,m) - obj.vertex(:,vop_p) ) , ( obj.vertex(:,vp_2) - obj.vertex(:,vop_p) ) );
    obj_V.un_ad( : , (m-1)*5 + 2 )  = unitary( term_cross_2 );
    obj_V.ds_ad( (m-1)*5 + 2 ) = sqrt( sum(term_cross_2.^2) )./2;
    obj_V.cent_ad( : , (m-1)*5 + 2 ) = (1/3)*( obj.vertex(:,vop_p) + obj.vertex(:,vp_2) + obj_V.vert_ad(:,m) );    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% [Third] Central triangle %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    obj_V.topol_ad( : , (m-1)*5 + 3 ) = [ vp_1   vp_2   (-m) ];
    term_cross_3  = cross( ( obj.vertex(:,vp_2) - obj.vertex(:,vp_1) ) , ( obj_V.vert_ad(:,m) - obj.vertex(:,vp_1) ) ); 
    obj_V.un_ad( : , (m-1)*5 + 3 )  = unitary( term_cross_3 );
    obj_V.ds_ad( (m-1)*5 + 3 ) = sqrt( sum(term_cross_3.^2) )./2;
    obj_V.ds_ad_nrm( m ) = (1/2) * obj.ln(m); 
    obj_V.cent_ad(: , (m-1)*5 + 3 ) = (1/3)*( obj.vertex(:,vp_1) + obj.vertex(:,vp_2) + obj_V.vert_ad(:,m) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% [Fourth] Added-negative triangle %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    obj_V.topol_ad( : , (m-1)*5 + 4 ) = [ vm_1   vop_m  (-m) ];
    term_cross_4 = cross( ( obj.vertex(:,vop_m) - obj.vertex(:,vm_1) ) , ( obj_V.vert_ad(:,m) - obj.vertex(:,vm_1) ) );    
    obj_V.un_ad( : , (m-1)*5 + 4 )  = unitary( term_cross_4 );
    obj_V.ds_ad( (m-1)*5 + 4 ) = sqrt( sum(term_cross_4.^2) )./2;
    obj_V.cent_ad( : , (m-1)*5 + 4 ) = (1/3)*( obj.vertex(:,vm_1) + obj.vertex(:,vop_m) + obj_V.vert_ad(:,m) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%% [Fifth] Added-negative triangle %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obj_V.topol_ad( : , (m-1)*5 + 5 ) = [ vop_m   vm_2   (-m) ];
    term_cross_5 = cross( ( obj.vertex(:,vm_2) - obj.vertex(:,vop_m) ) , ( obj_V.vert_ad(:,m) - obj.vertex(:,vop_m) ) );        
    obj_V.un_ad( : , (m-1)*5 + 5 )  = unitary( term_cross_5 );
    obj_V.ds_ad( (m-1)*5 + 5 ) = sqrt( sum(term_cross_5.^2) )./2;
    obj_V.cent_ad( : , (m-1)*5 + 5 ) = (1/3)*( obj.vertex(:,vm_2) + obj.vertex(:,vop_m) + obj_V.vert_ad(:,m) );
    
end;  %% for m=1:Ne,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    (II) VOLUME PARAMETERS    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_V.monoSWG = zeros( 8 , 2*Ne );
obj_V.cent_V  = zeros( 3 , 2*Ne );
obj_V.dV  = zeros( 1 , 2*Ne );
obj_V.dV_nrm  = zeros( 1 , 2*Ne );

for m=1:Ne,
    
    %%%% Volume plus update
    Vp = m;
    Tp = obj.edges(1,m);
    vp_or = obj.topol(:,Tp);
    
    %%% obj_V Triangles update
    obj_V.monoSWG( 1:4 , Vp) = [ Tp  (-((m-1)*5 + [1:3])) ].';
    %%% obj_V Vertices update
    obj_V.monoSWG( 5:8 , Vp) = [ vp_or.'  (-m)].';
    
    %%%% obj_V Volume minus update
    Vm = Ne + m;
    Tm = obj.edges(2,m);
    vm_or = obj.topol(:,Tm);    
    
    %%% obj_V Triangles update
    obj_V.monoSWG( 1:4 , Vm ) = [ Tm  (-((m-1)*5 + [3:5])) ].';
    %%% Vertices update
    obj_V.monoSWG( 5:8 , Vm ) = [ vm_or.' (-m)].';
    
    %%% update central triangle
    obj_V.Tmid( m ) = (m-1)*5 + 3; 
    %
    obj_V.Tp_ad1( m ) = (m-1)*5 + 1; 
    obj_V.Tp_ad2( m ) = (m-1)*5 + 2; 
    %
    obj_V.Tm_ad1( m ) = (m-1)*5 + 4;     
    obj_V.Tm_ad2( m ) = (m-1)*5 + 5;         
    
    %%%% Obj_V centroid plus
    obj_V.cent_V( : , Vp ) = ( sum( obj.vertex(:,vp_or) .' ) .' + obj_V.vert_ad(:,m) )./4; 
    %%%% Obj_V centroid minus
    obj_V.cent_V( : , Vm ) = ( sum( obj.vertex(:,vm_or) .' ) .' + obj_V.vert_ad(:,m) )./4;     
    
    %%%% obj_V volume plus
    obj_V.dV( Vp ) = (1/3) * obj.ds(Tp) * depth_length(m) * sum( ( -obj.un(:,Tp) ) .* obj_V.un_ed_int(:,m) );
    obj_V.dV_nrm( Vp ) = (1/3) * obj.ds(Tp) * sum( ( -obj.un(:,Tp) ) .* obj_V.un_ed_int(:,m) );
    
    %%%% obj_V volume minus
    obj_V.dV( Vm ) = (1/3) * obj.ds(Tm) * depth_length(m) * sum( ( -obj.un(:,Tm) ) .* obj_V.un_ed_int(:,m) );
    obj_V.dV_nrm( Vm ) = (1/3) * obj.ds(Tm) * sum( ( -obj.un(:,Tm) ) .* obj_V.un_ed_int(:,m) );    
    
end;  %%% for m=1:Ne,
