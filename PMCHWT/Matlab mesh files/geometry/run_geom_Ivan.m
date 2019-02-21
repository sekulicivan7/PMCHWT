%  Generates geometry
%
%   field==1       --> EFIE[R] [ed]
%   field==2       --> MFIE[R] [ed]
%   field==3       --> MFIE[monoR] [ed]
%   field==4       --> EFIE[tet][monoR] [tri]
%   field==5       --> EFIE[tet][monoR] [ed]
%   field==6       --> EFIE[wed][monoR] [tri]
%   field==7       --> EFIE[wed][monoR] [ed]
%   field==8       --> EFIE[tet][eoR] [ed]
%   field==9       --> EFIE[wed][eoR] [ed]
%   field==10       --> EFIE[TN][monoR] [tri]
%   field==11       --> EFIE[TN][monoR] [ed]
%
% by Ed. Ubeda [February 2016]

disp('Processing geometry:');
tic;

if (param==0),
    %%% gid object file
    p_obj.gid.file = geom;
    m_gidmesh_Ivan;
else,
    %%% parameters
    
    obj= feval(geom,param);
    obj = get_edge_Ivan(obj);
end;

%%% Gram matrices
%% D = <RWG,RWG>
D = d_unnrm_mat(obj);

%% D_mono = <monoR,monoR> (edge-reordered)
D_mono = d_monorwg_mat(obj);

%%% D for even-odd MFIE[tet/wed]
[D_e_e, D_e_o] = d_unnrm_ee_eo_mat(obj);

%% D_tri_mono = <monoR,monoR> (triangle-reordered)
Ne = size(obj.edges,2);
Nt = size(obj.topol,2);

set_ed_v1 = obj.trian(1,:);
i_neg_v1 = find(set_ed_v1<0);
set_ed_v1(i_neg_v1) = abs(set_ed_v1(i_neg_v1)) + Ne;

set_ed_v2 = obj.trian(2,:);
i_neg_v2 = find(set_ed_v2<0);
set_ed_v2(i_neg_v2) = abs(set_ed_v2(i_neg_v2)) + Ne;

set_ed_v3 = obj.trian(3,:);
i_neg_v3 = find(set_ed_v3<0);
set_ed_v3(i_neg_v3) = abs(set_ed_v3(i_neg_v3)) + Ne;

%%% set_ed_2_tri
set_ed_2_tri = zeros(1,3*Nt);
set_ed_2_tri(1:3:end) = set_ed_v1;
set_ed_2_tri(2:3:end) = set_ed_v2;
set_ed_2_tri(3:3:end) = set_ed_v3;        

D_tri_mono = D_mono(set_ed_2_tri,set_ed_2_tri);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TETRAHEDRAL GEOMETRIC PARAMETERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ l_vec, vert_or, vert_mid, un_c_plus, un_c_min ] = get_l_vec(obj);
obj_V = get_obj_vol_Ivan( obj , un_c_plus , un_c_min , vert_mid , l_vec , fact_vol );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOP GEOMETRIC PARAMETERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('obj'),

    disp('Computing Loop: Original mesh');
    tic;

    loop_or = struct('tri',[],'ed',[],'edc',[],'edc_0',[],'vc',[],'vc_0',[],'nv',[],'ntriv',[],'tri2ver',[],'tri_drain',[],'ver2tri',[],'tri2ver_2',[]);

    %%% Find the number maximum of items around a vertex arising in the discretization
    Nv_or = length( obj.vertex );
    Nt_or = length( obj.ds );

    n_acc_tmp_or = 0;

    for p=1:Nv_or,
        tmp = find( (p==obj.topol(1,:)) | (p==obj.topol(2,:)) | (p==obj.topol(3,:)) );
        n_tmp = length(tmp);

        n_acc_tmp_or = max([n_acc_tmp_or n_tmp]);
    end;    %%% for p=1:Nv_or,

    %%%% Setup of loop_or
    loop_or.tri = zeros(Nv_or,n_acc_tmp_or);
    loop_or.ed = zeros(Nv_or,n_acc_tmp_or);
    loop_or.edc = zeros(Nv_or,n_acc_tmp_or);
    loop_or.edc_0 = zeros(Nv_or,n_acc_tmp_or);    
    loop_or.vc = zeros(Nv_or,n_acc_tmp_or);
    loop_or.vc_0 = zeros(Nv_or,n_acc_tmp_or);
    loop_or.tri2ver = sparse(Nt_or,Nv_or);
    loop_or.tri2ver_2 = sparse(Nt_or,Nv_or);
    loop_or.ver2tri = zeros(Nv_or,n_acc_tmp_or);
    %%% accumulated area of the triangles around the vertex
    loop_or.ds = zeros(1,Nv_or);

    loop_or.nv = zeros(Nv_or,1);

    loop_or.ntriv = [];

    %%%% Computation of loop_or:
    for p=1:Nv_or,

        %% update loop_or.tri
        tmp_T = find( (p==obj.topol(1,:)) | (p==obj.topol(2,:)) | (p==obj.topol(3,:)) );
        tmp_end_T = [];

        tmp_end_vc = [];
        tmp_end_vc_alt = [];

        %% Right ordering of tmp_T
        done=0;

        %%% First triangle
        T0 = tmp_T(1);

        v0 = obj.topol(:,T0);
        i_v_0 = find(v0~=p);

        i_v_yes = find(v0==p);

        if ( ( (i_v_0(1))==1 & (i_v_0(2)==2) ) | ( (i_v_0(2))==1 & (i_v_0(1)==2) ) ),
            v_0 = v0(2);
            vc_0 = v0(1);
        end; %%% if ( ( (i_v_0(1))==1 & (i_v_0(2)==2) ) | ( (i_v_0(2))==1 & (i_v_0(1)==2) ) ),
        %
        if ( ( (i_v_0(1))==1 & (i_v_0(2)==3) ) | ( (i_v_0(2))==1 & (i_v_0(1)==3) ) ),
            v_0 = v0(1);
            vc_0 = v0(3);
        end; %%% if ( ( (i_v_0(1))==1 & (i_v_0(2)==3) ) | ( (i_v_0(2))==1 & (i_v_0(1)==3) ) ),
        %
        if ( ( (i_v_0(1))==2 & (i_v_0(2)==3) ) | ( (i_v_0(2))==2 & (i_v_0(1)==3) ) ),
            v_0 = v0(3);
            vc_0 = v0(2);
        end; %%% if ( ( (i_v_0(1))==1 & (i_v_0(2)==3) ) | ( (i_v_0(2))==1 & (i_v_0(1)==3) ) ),

        tmp_end_T = T0;
        tmp_end_vc = v_0;

        %% new
        tmp_end_vc_alt = vc_0;
        tmp_end_v_yes = i_v_yes;
        %%

        tmp_T(1)=0;

        while (done==0),

            for q=1:length(tmp_T),

                if (tmp_T(q)~=0),

                    v_q = obj.topol(:,tmp_T(q));
                    i_v_q = find(v_q~=p);
                    v_q = v_q(i_v_q);

                    %% Comparison
                    i_v_q = find( v_q == v_0 );
                    if length(i_v_q),
                        i_v_noq = find( v_q ~= v_0 );
                        vc_0  = v_0;
                        v_0 = v_q(i_v_noq);
                        tmp_end_T = [ tmp_end_T tmp_T(q) ];
                        tmp_end_vc = [ tmp_end_vc v_0 ];
                        %% new
                        tmp_end_vc_alt = [ tmp_end_vc_alt vc_0];
                        %%
                        %% find local vertex index in the triangle
                        i_v_yes = find(p==obj.topol(:,tmp_T(q)));
                        tmp_end_v_yes = [ tmp_end_v_yes   i_v_yes ];
                        %%
                        tmp_T(q) = 0;
                        break;
                    end;  %%% if length(i_v_q),

                end;  %%% if (tmp_T~=0),

            end; %%% for q=1:(length(tmp_T)-i0),

            if ( length(tmp_end_T)==length(tmp_T) ), done=1;  end;

        end;  %% while (done==0),

        loop_or.tri(p,1:length(tmp_T)) = tmp_end_T;
        loop_or.ds(p) = sum( obj.ds(tmp_end_T).' );

        loop_or.ver2tri(p,1:length(tmp_end_v_yes)) = tmp_end_v_yes;
        loop_or.tri_drain(p) = tmp_end_T(end);
        loop_or.ver2tri_drain(p) = tmp_end_v_yes(end);
        loop_or.tri2ver(tmp_end_T,p) = (1:length(tmp_end_T)).';
        if p~=Nv_or,
            %% the Nv_or vertex is excluded in the loop basis functions
            loop_or.tri2ver_2(tmp_end_T,p) = (1:length(tmp_end_T)).';
        end;
        %% TAG the drain-triangle: the last one
        %%%% IF REQUIRED THE ORDER GETS FLIPPED
        %%%% loop_or.tri2ver(tmp_end_T(end),p) = - loop_or.tri2ver(tmp_end_T(end),p);
        %%
        %%
        loop_or.vc(p,1:length(tmp_T)) = tmp_end_vc;
        %% new
        loop_or.vc_0(p,1:length(tmp_T)) = tmp_end_vc_alt;
        
        %%% Find adjacent edges
        tmp_end_edc = [];
        for q=1:length(tmp_end_T),
            i_edc = find( ( ( tmp_end_T(q)==obj.edges(1,:) ) & ( tmp_end_vc(q)==obj.edges(3,:) ) ) | ...
                ( ( tmp_end_T(q)==obj.edges(2,:) ) & ( tmp_end_vc(q)==obj.edges(4,:) ) ) );
            tmp_end_edc = [tmp_end_edc i_edc];
        end;    %% for q=1:length(tmp_end_T),

        loop_or.edc( p , 1:length(tmp_T) ) = tmp_end_edc;
        loop_or.edc_0( p , 1:length(tmp_T) ) = [ tmp_end_edc(end) tmp_end_edc(1:end-1) ];
        
        %%% Find oposed edges
        tmp_end_ed = [];
        for q=1:length(tmp_end_T),
            i_ed = find( ( ( tmp_end_T(q)==obj.edges(1,:) ) & ( p==obj.edges(3,:) ) ) | ...
                ( ( tmp_end_T(q)==obj.edges(2,:) ) & ( p==obj.edges(4,:) ) ) );
            tmp_end_ed = [tmp_end_ed i_ed];
        end;    %% for q=1:length(tmp_end_T),
        
        loop_or.ed(p,1:length(tmp_T)) = tmp_end_ed;

        %%% Number of items per vertex
        %%% ===>
        %%% ===>
        l_T = length(tmp_T);
        loop_or.nv(p) = l_T;
        %%% <===
        %%% <===

        %%% Number of triangles around vertices
        if (length(loop_or.ntriv)==0),
            loop_or.ntriv = l_T;
        else,
            if all(l_T~=loop_or.ntriv),
                loop_or.ntriv = [loop_or.ntriv l_T];
            end;    %% if length( find(l_T~=loop_or.ntriv) ),
        end;    %% if (~length(loop_or.ntriv)),

    end; %%% for p=1:Nv_or,
    
    %%% Obtain de Star basis functions
    [star_or, tri_star_or] = get_star_d3(obj);

    disp(sprintf('time = %g',toc));

end;  %%% if exist('obj'),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WEDGE GEOMETRIC PARAMETERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Processing INSIDE geometry:');
obj_inpV = get_obj_inprismV_Ivan( obj , obj_V, loop_or , fact_vol );
obj_nc_inpV = get_obj_nc_inprismV_Ivan( obj , fact_vol );
obj_nc_inpV2 = get_obj_nc_inprismV_Ivan2( obj , fact_vol );
obj_inpV2 = get_obj_inprismV_Ivan2( obj , obj_V, loop_or , fact_vol );

disp(sprintf('time = %g',toc));