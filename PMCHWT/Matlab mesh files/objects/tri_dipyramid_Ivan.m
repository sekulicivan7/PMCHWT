% A regular triangular dipyramid is built
%
% param(1)  -----> side length
% param(2)  -----> steps on the side
%
% The symmetric base of the triangular dipyramid lies on the plane XY
% [This is a copied tretrahedron with simmetry wrt the XY plane]
%
% by Eduard Ubeda Farre, September 2004.
%

function obj = tri_dipyramid_Ivan(param)

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'edges_c',[],'un',[],'ds',[],'ln',[],'cent',[]);

L = param(1);
N = param(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the structure of the TETRAHEDRON %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of triangles per face
nt = sum( [1 : 2 : (2*N-1)].' );

obj.topol = zeros( 3 , 4*nt );

%%%% vertex_xpl_1= [ L/sqrt(3) 0 0 ].';
vertex_xpl_1= [0 0 sqrt(2/3)*L].';
vertex_xmin_ymin_2 = [ -L/(2*sqrt(3))  -L/2  0 ].';
vertex_xmin_ypl_3  = [ -L/(2*sqrt(3))  +L/2  0 ].';
%%%% vertex_zpl_4 = [0 0 sqrt(2/3)*L].';
vertex_zpl_4 = [ L/sqrt(3) 0 0 ].';

%%% new for hexa_hed_d
%%% ==>
%%%%% vertex_xmin_5= [ -L/sqrt(3) 0 0 ].';
vertex_xmin_5= [0 0 -sqrt(2/3)*L].';
%%% <==
vc_1_2 = vertex_xmin_ymin_2 - vertex_xpl_1;
l_1_2 = norm( vc_1_2 );

vc_1_3 = vertex_xmin_ypl_3 - vertex_xpl_1;
l_1_3 = norm( vc_1_3 );

vc_1_4 = vertex_zpl_4 - vertex_xpl_1;
l_1_4 = norm( vc_1_4 );

vc_2_3 = vertex_xmin_ypl_3 - vertex_xmin_ymin_2;
l_2_3 = norm( vc_2_3 );

vc_2_4 = vertex_zpl_4 - vertex_xmin_ymin_2;
l_2_4 = norm( vc_2_4 );

vc_3_4 = vertex_zpl_4 - vertex_xmin_ypl_3;
l_3_4 = norm( vc_3_4 );

%%% new for hexa_hed_d
%%% ==>
vc_5_2 = vertex_xmin_ymin_2 - vertex_xmin_5;
l_5_2 = norm( vc_1_2 );

vc_5_3 = vertex_xmin_ypl_3 - vertex_xmin_5;
l_5_3 = norm( vc_5_3 );

vc_5_4 = vertex_zpl_4 - vertex_xmin_5;
l_5_4 = norm( vc_5_4 );
%%% <==

set = 0:(1/N):1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINITION OF VERTICES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (I)  Define vertices ON VERTICES. Update

%%%%%%%%%%%%%%%%%%%%
%%%% UPPER HALF %%%%
%%%%%%%%%%%%%%%%%%%%
%%%% ==>>>
%%%% ==>>>
obj.vertex = [ vertex_xpl_1  vertex_xmin_ymin_2  vertex_xmin_ypl_3  vertex_zpl_4 ];

fst_vert_1_2_3 = 1;
fst_vert_1_3_4 = 1;
fst_vert_1_4_2 = 1;
fst_vert_2_3_4 = 2;

n_vert_2 = 2;
n_vert_3 = 3;
n_vert_4 = 4;

%% (II) Define vertices ON EDGES
vert_edge_1_2 = vertex_xpl_1*ones(size(set))  +  vc_1_2 * set;
vert_edge_1_3 = vertex_xpl_1*ones(size(set))  +  vc_1_3 * set;
vert_edge_1_4 = vertex_xpl_1*ones(size(set))  +  vc_1_4 * set;

vert_edge_2_3 = vertex_xmin_ymin_2*ones(size(set))  +  vc_2_3 * set;
vert_edge_2_4 = vertex_xmin_ymin_2*ones(size(set))  +  vc_2_4 * set;

vert_edge_3_4 = vertex_xmin_ypl_3*ones(size(set))  +  vc_3_4 * set;

%% update
obj.vertex = [ obj.vertex vert_edge_1_2(:,2:N)  vert_edge_1_3(:,2:N)  vert_edge_1_4(:,2:N) ...
                          vert_edge_2_3(:,2:N)  vert_edge_2_4(:,2:N)  vert_edge_3_4(:,2:N)  ];
                  
vert_line_1_2 = 4 + (1:(N-1));
vert_line_1_3 = 4 + (N-1) + (1:(N-1));
vert_line_1_4 = 4 + 2*(N-1) + (1:(N-1));
vert_line_2_3 = 4 + 3*(N-1) + (1:(N-1));
vert_line_2_4 = 4 + 4*(N-1) + (1:(N-1));
vert_line_4_2 = fliplr( vert_line_2_4 );
vert_line_3_4 = 4 + 5*(N-1) + (1:(N-1));

%% (III) Define vertices INSIDE FACES
%% The convention for filling the faces x_y_z of the triangles is from edge x_y to edge x_z 
%% ( parallel to edge x ) One must ensure that the vectorial product is normal to the surface 
%% with the same sense. Define vertices on faces
%% 1st face : 1_2_3
%% 2nd face : 1_3_4
%% 3rd face : 1_4_2
%% 4th face : 2_3_4
%% update

n_1_2_3_0 = length( obj.vertex );

vert_face_1_2_3 = [];
for m=3:N,
    tmp_vec = ( vert_edge_1_3(:,m) - vert_edge_1_2(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_1_2(:,m)*ones(1,(m-2)) + tmp_vec * [ 1 : (m-2) ];
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

n_1_3_4_0 = length( obj.vertex );

vert_face_1_3_4 = [];
for m=3:N,
    tmp_vec = ( vert_edge_1_4(:,m) - vert_edge_1_3(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_1_3(:,m)*ones(1,(m-2)) + tmp_vec * [ 1 : (m-2) ];
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

n_1_4_2_0 = length( obj.vertex );

vert_face_1_4_2 = [];
for m=3:N,
    tmp_vec = ( vert_edge_1_2(:,m) - vert_edge_1_4(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_1_4(:,m)*ones(1,m-2) +   tmp_vec *[ 1 : (m-2) ] ;
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

%%%%%% <====
%%%%%% <====

%%%%%%%%%%%%%%%%%%%%
%%%% LOWER HALF %%%%
%%%%%%%%%%%%%%%%%%%%
%%%% new for hexa_hed_d
%%%% ==>>>
%%%% ==>>>

Nv_up = length( obj.vertex );

%%% remaining vertices:  vertex_xmin_ymin_2  vertex_xmin_ypl_3 vertex_zpl_4
obj.vertex = [ obj.vertex vertex_xmin_5  ];

fst_vert_5_2_3 = Nv_up + 1;
fst_vert_5_3_4 = Nv_up + 1;
fst_vert_5_4_2 = Nv_up + 1;

%%%% n_vert_2 = 2;
%%%% n_vert_3 = 3;
%%%% n_vert_4 = 4;
n_vert_5 = Nv_up + 1;

%% (II) Define vertices ON EDGES
vert_edge_5_2 = vertex_xmin_5*ones(size(set))  +  vc_5_2 * set;
vert_edge_5_3 = vertex_xmin_5*ones(size(set))  +  vc_5_3 * set;
vert_edge_5_4 = vertex_xmin_5*ones(size(set))  +  vc_5_4 * set;

%%%% vert_edge_2_3 = vertex_xmin_ymin_2*ones(size(set))  +  vc_2_3 * set;
%%%% vert_edge_2_4 = vertex_xmin_ymin_2*ones(size(set))  +  vc_2_4 * set;

%%%% vert_edge_3_4 = vertex_xmin_ypl_3*ones(size(set))  +  vc_3_4 * set;

%% update
obj.vertex = [ obj.vertex vert_edge_5_2(:,2:N)  vert_edge_5_3(:,2:N)  vert_edge_5_4(:,2:N) ];

vert_line_5_2 = Nv_up + 1 + (1:(N-1));
vert_line_5_3 = Nv_up + 1 + (N-1) + (1:(N-1));
vert_line_5_4 = Nv_up + 1 + 2*(N-1) + (1:(N-1));
%%%% vert_line_2_3 = 4 + 3*(N-1) + (1:(N-1));
%%%% vert_line_2_4 = 4 + 4*(N-1) + (1:(N-1));
%%%% vert_line_4_2 = fliplr( vert_line_2_4 );
%%%% vert_line_3_4 = 4 + 5*(N-1) + (1:(N-1));

%% (III) Define vertices INSIDE FACES
%% The convention for filling the faces x_y_z of the triangles is from edge x_y to edge x_z 
%% ( parallel to edge x ) One must ensure that the vectorial product is normal to the surface 
%% with the same sense. Define vertices on faces
%% 1st face : 1_2_3
%% 2nd face : 1_3_4
%% 3rd face : 1_4_2
%% 4th face : 2_3_4
%% update

n_5_2_3_0 = length( obj.vertex );

vert_face_5_2_3 = [];
for m=3:N,
    tmp_vec = ( vert_edge_5_3(:,m) - vert_edge_5_2(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_5_2(:,m)*ones(1,(m-2)) + tmp_vec * [ 1 : (m-2) ];
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

n_5_3_4_0 = length( obj.vertex );

vert_face_5_3_4 = [];
for m=3:N,
    tmp_vec = ( vert_edge_5_4(:,m) - vert_edge_5_3(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_5_3(:,m)*ones(1,(m-2)) + tmp_vec * [ 1 : (m-2) ];
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

n_5_4_2_0 = length( obj.vertex );

vert_face_5_4_2 = [];
for m=3:N,
    tmp_vec = ( vert_edge_5_2(:,m) - vert_edge_5_4(:,m) )/(m-1); 
    vert_edge_tmp = vert_edge_5_4(:,m)*ones(1,m-2) +   tmp_vec *[ 1 : (m-2) ] ;
    obj.vertex = [ obj.vertex  vert_edge_tmp ];
end;

%%%%%% <====
%%%%%% <====

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DEFINITION OF TRIANGLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% UPPER HALF %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ===>
%%%%% ===>
obj.topol = [];
n_prev_1 = n_1_2_3_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 1st face: 1_2_3  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st row 
tmp_topol = [ fst_vert_1_2_3  vert_line_1_2(1)  vert_line_1_3(1)   ].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_1_2(1)   vert_line_1_2(2) (n_prev_1 + 1) ].' ...
                        [ vert_line_1_3(1)   vert_line_1_2(1)  (n_prev_1 + 1)     ].' ...
                        [ vert_line_1_3(1)   (n_prev_1 + 1) vert_line_1_3(2)    ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_1_2(m-1)   vert_line_1_2(m) (n_prev_1 + 1) ].' ...
                                [ vert_line_1_2(m-1)   (n_prev_1 + 1)  (n_prev_0 + 1)  ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n)  (n_prev_1 + n)   (n_prev_1 + (n+1))  ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_0 + n) (n_prev_1 + n+1) (n_prev_0 + n+1)  ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [       n_prev_1    (n_prev_1 + (m-1))  vert_line_1_3(m-1)  ].' ...
                                [ (n_prev_1 + (m-1))   vert_line_1_3(m)  vert_line_1_3(m-1)    ].'  ];                
                        
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_2_3, n_vert_2, n_vert_3) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
                                       
        tmp_topol = [ tmp_topol [ vert_line_1_2(N-1) n_vert_2  vert_line_2_3(1)  ].' ...
                                [ vert_line_1_2(N-1) vert_line_2_3(1)  (n_prev_0+1)  ].'     ];
                       
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [ vert_line_2_3(n) vert_line_2_3(n+1) (n_prev_0 + n)  ].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  (n_prev_0 + n) vert_line_2_3(n+1) (n_prev_0 + (n+1)) ].'  ]; 
            end;
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  vert_line_2_3(N-1)  vert_line_1_3(N-1) (n_prev_0 + (N-2)) ].'     ...
                                    [  vert_line_1_3(N-1)  vert_line_2_3(N-1)    n_vert_3        ].'    ];
                            
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

n_prev_1 = n_1_3_4_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 2nd face: 1_3_4  %% 1_2_3 => 1_3_4 ; 1_2 => 1_3 ; 1_3 => 1_4 ; 2_3 => 3_4 ; n_vert_2 => n_vert_3 ; n_vert_3 => n_vert_4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st row 
tmp_topol = [ fst_vert_1_3_4  vert_line_1_3(1)  vert_line_1_4(1) ].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_1_3(1)  vert_line_1_3(2) (n_prev_1 + 1) ].' ...
                        [ vert_line_1_3(1)  (n_prev_1 + 1)  vert_line_1_4(1)  ].' ...
                        [ vert_line_1_4(1)   (n_prev_1 + 1) vert_line_1_4(2) ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_1_3(m-1)  vert_line_1_3(m) (n_prev_1 + 1) ].' ...
                                [ vert_line_1_3(m-1) (n_prev_1 + 1)  (n_prev_0 + 1) ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n) (n_prev_1 + n)  (n_prev_1 + (n+1))  ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_0 + n) (n_prev_1 + n+1) (n_prev_0 + n+1) ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [ (n_prev_1 + (m-1))    vert_line_1_4(m-1)    n_prev_1 ].' ...
                                [ (n_prev_1 + (m-1))    vert_line_1_4(m)   vert_line_1_4(m-1) ].'  ];                
                        
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_3_4, n_vert_3, n_vert_4) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
                                       
        tmp_topol = [ tmp_topol [ vert_line_1_3(N-1)  n_vert_3  vert_line_3_4(1) ].' ...
                                [ vert_line_1_3(N-1)   vert_line_3_4(1) (n_prev_0+1) ].'     ];
                       
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [  vert_line_3_4(n) vert_line_3_4(n+1) (n_prev_0 + n) ].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  vert_line_3_4(n+1)  (n_prev_0 + (n+1)) (n_prev_0 + n) ].'  ]; 
            end; 
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  (n_prev_0 + (N-2))  vert_line_3_4(N-1)   vert_line_1_4(N-1)  ].'     ...
                                    [  vert_line_1_4(N-1)  vert_line_3_4(N-1)    n_vert_4         ].'    ];
                            
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

n_prev_1 = n_1_4_2_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 3rd face: 1_4_2  %% 1_2_3 => 1_4_2 ;  1_2 => 1_4 ; 1_3 => 1_2 ; 2_3 => 4_2 ; n_vert_2 => n_vert_4 ; n_vert_3 => n_vert_2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st row 
tmp_topol = [ fst_vert_1_4_2  vert_line_1_4(1) vert_line_1_2(1)  ].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_1_4(1)  vert_line_1_4(2)  (n_prev_1 + 1)  ].' ...
                        [ vert_line_1_4(1)  (n_prev_1 + 1)  vert_line_1_2(1)   ].' ...
                        [ vert_line_1_2(1)  (n_prev_1 + 1) vert_line_1_2(2)  ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_1_4(m)  (n_prev_1 + 1)  vert_line_1_4(m-1)  ].' ...
                                [  (n_prev_1 + 1)   (n_prev_0 + 1)  vert_line_1_4(m-1)  ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_1 + n)  (n_prev_1 + (n+1))  (n_prev_0 + n)  ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_1 + n+1)  (n_prev_0 + n+1) (n_prev_0 + n)  ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [   (n_prev_1 + (m-1))   vert_line_1_2(m-1)     n_prev_1          ].' ...
                                [   vert_line_1_2(m)   vert_line_1_2(m-1)     (n_prev_1 + (m-1))   ].'  ];                
                        
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_4_2, n_vert_4, n_vert_2) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
                                       
        tmp_topol = [ tmp_topol [ vert_line_1_4(N-1)  n_vert_4  vert_line_4_2(1) ].' ...
                                [ vert_line_1_4(N-1)  vert_line_4_2(1)  (n_prev_0+1)    ].'     ];
                       
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n)  vert_line_4_2(n)  vert_line_4_2(n+1) ].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  vert_line_4_2(n+1)  (n_prev_0 + (n+1))  (n_prev_0 + n)  ].'  ]; 
            end;
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  (n_prev_0 + (N-2))  vert_line_4_2(N-1)  vert_line_1_2(N-1)   ].'     ...
                                    [  vert_line_1_2(N-1)  vert_line_4_2(N-1)     n_vert_2         ].'    ];
                            
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

%%%%% <===
%%%%% <===

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOWER HALF %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ===>
%%%%% ===>
n_prev_1 = n_5_2_3_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 4th face: 5_2_3  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st row 
tmp_topol = [ fst_vert_5_2_3    vert_line_5_3(1)  vert_line_5_2(1) ].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_5_2(1)   (n_prev_1 + 1)   vert_line_5_2(2) ].' ...
                        [ vert_line_5_3(1)   (n_prev_1 + 1)   vert_line_5_2(1)  ].' ...
                        [ vert_line_5_3(1)   vert_line_5_3(2)   (n_prev_1 + 1) ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_5_2(m-1)   (n_prev_1 + 1)  vert_line_5_2(m)].' ...
                                [ vert_line_5_2(m-1)   (n_prev_0 + 1)  (n_prev_1 + 1)    ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n)  (n_prev_1 + (n+1))  (n_prev_1 + n) ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_0 + n) (n_prev_0 + n+1) (n_prev_1 + n+1) ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [       n_prev_1     vert_line_5_3(m-1) (n_prev_1 + (m-1)) ].' ...
                                [ (n_prev_1 + (m-1))  vert_line_5_3(m-1)  vert_line_5_3(m)  ].'  ];                
                            
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_2_3, n_vert_2, n_vert_3) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
        tmp_topol = [ tmp_topol [ vert_line_5_2(N-1)   vert_line_2_3(1) n_vert_2 ].' ...
                                [ vert_line_5_2(N-1) (n_prev_0+1) vert_line_2_3(1) ].'     ];
                            
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [ vert_line_2_3(n) (n_prev_0 + n)   vert_line_2_3(n+1) ].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  (n_prev_0 + n)  (n_prev_0 + (n+1))  vert_line_2_3(n+1) ].'  ]; 
            end;
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  vert_line_2_3(N-1)  (n_prev_0 + (N-2))  vert_line_5_3(N-1)   ].'     ...
                                    [  vert_line_5_3(N-1)       n_vert_3       vert_line_2_3(N-1)   ].'    ];
                                
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

n_prev_1 = n_5_3_4_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 5th face: 5_3_4  %% 5_2_3 => 5_3_4 ; 5_2 => 5_3 ; 5_3 => 5_4 ; 5_3 => 5_4 ; n_vert_2 => n_vert_3 ; n_vert_3 => n_vert_4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st row 
tmp_topol = [ fst_vert_5_3_4  vert_line_5_4(1)  vert_line_5_3(1)].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_5_3(1)  (n_prev_1 + 1)  vert_line_5_3(2) ].' ...
                        [ vert_line_5_3(1)  vert_line_5_4(1)  (n_prev_1 + 1) ].' ...
                        [ vert_line_5_4(1)  vert_line_5_4(2)  (n_prev_1 + 1) ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_5_3(m-1)  (n_prev_1 + 1)  vert_line_5_3(m) ].' ...
                                [ vert_line_5_3(m-1)  (n_prev_0 + 1) (n_prev_1 + 1)  ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n) (n_prev_1 + (n+1)) (n_prev_1 + n)  ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_0 + n)   (n_prev_0 + n+1) (n_prev_1 + n+1) ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [ (n_prev_1 + (m-1))        n_prev_1        vert_line_5_4(m-1) ].' ...
                                [ (n_prev_1 + (m-1))    vert_line_5_4(m-1)   vert_line_5_4(m)   ].'  ];                
                        
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_3_4, n_vert_3, n_vert_4) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
        tmp_topol = [ tmp_topol [ vert_line_5_3(N-1)  vert_line_3_4(1) n_vert_3].' ...
                    [ vert_line_5_3(N-1)   (n_prev_0+1)  vert_line_3_4(1)].'     ];
        
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [  vert_line_3_4(n) (n_prev_0 + n)   vert_line_3_4(n+1)].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  vert_line_3_4(n+1)  (n_prev_0 + n)  (n_prev_0 + (n+1)) ].'  ]; 
            end; 
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  (n_prev_0 + (N-2))    vert_line_5_4(N-1) vert_line_3_4(N-1)  ].'     ...
                                    [  vert_line_5_4(N-1)      n_vert_4        vert_line_3_4(N-1)   ].'    ];
                            
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

n_prev_1 = n_5_4_2_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  TRIANGLE : 6th face: 5_4_2  %% 5_2_3 => 5_4_2 ;  5_2 => 5_4 ; 5_3 => 5_2 ; 5_3 => 4_2 ; n_vert_2 => n_vert_4 ; n_vert_3 => n_vert_2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st row 
tmp_topol = [ fst_vert_5_4_2  vert_line_5_2(1) vert_line_5_4(1)  ].';
% 2nd row
tmp_topol = [ tmp_topol [ vert_line_5_4(1)   (n_prev_1 + 1)   vert_line_5_4(2) ].' ...
                        [ vert_line_5_4(1)  vert_line_5_2(1)  (n_prev_1 + 1) ].' ...
                        [ vert_line_5_2(1)  vert_line_5_2(2)  (n_prev_1 + 1) ].' ];
               
%% Remanining  N-2 rows
n_prev_0 = n_prev_1;
n_prev_1 = n_prev_0 + 1;

for m=3:N
    
    if (m~=N),
        
        %% First-row triangles
        tmp_topol = [ tmp_topol [ vert_line_5_4(m)   vert_line_5_4(m-1)  (n_prev_1 + 1) ].' ...
                                [  (n_prev_1 + 1)    vert_line_5_4(m-1)  (n_prev_0 + 1) ].'   ];

        %% Mid-row triangles
        for n=1:m-2,
            tmp_topol = [ tmp_topol [ (n_prev_1 + n)  (n_prev_0 + n)  (n_prev_1 + (n+1)) ].' ];
            
            if ( n~=(m-2) ),                
                tmp_topol = [ tmp_topol [ (n_prev_1 + n+1)  (n_prev_0 + n)  (n_prev_0 + n+1) ].'  ];            
            end; %% if ( n~=(m-2) ),
            
        end;  %% for n=1:m-2,
                
        %% End-row triangles
        tmp_topol = [ tmp_topol [   (n_prev_1 + (m-1))     n_prev_1           vert_line_5_2(m-1) ].' ...
                                [   vert_line_5_2(m)   (n_prev_1 + (m-1))     vert_line_5_2(m-1) ].'  ];                
                        
        n_prev_0 = n_prev_1;
        n_prev_1 = n_prev_0 + (m-1);
        
    elseif (m==N),
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Last row of the triangle (vert_line_4_2, n_vert_4, n_vert_2) %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% First-row end-triangles
                                       
        tmp_topol = [ tmp_topol [ vert_line_5_4(N-1)  vert_line_4_2(1)  n_vert_4].' ...
                                [ vert_line_5_4(N-1)   (n_prev_0+1)   vert_line_4_2(1) ].'     ];
                       
        %% Mid-row end-triangles
        for n=1:N-2,
            tmp_topol = [ tmp_topol [ (n_prev_0 + n)  vert_line_4_2(n+1) vert_line_4_2(n) ].' ];
            
            if ( n~=(N-2) ),  
                tmp_topol = [ tmp_topol  [  vert_line_4_2(n+1)  (n_prev_0 + n)  (n_prev_0 + (n+1)) ].'  ]; 
            end;
        end;  %% for n=1:N-2,
        
        %% Last-row end-triangles
        tmp_topol = [ tmp_topol     [  (n_prev_0 + (N-2))  vert_line_5_2(N-1)  vert_line_4_2(N-1)   ].'     ...
                                    [  vert_line_5_2(N-1)      n_vert_2        vert_line_4_2(N-1)   ].'    ];
                            
    end;  %% if (m~=N),
    
end; %% for m=2:N

obj.topol= [ obj.topol tmp_topol ];
tmp_topol = [];

%%%%% <===
%%%%% <===
