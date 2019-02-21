% obj = sphere_Ivan(param)
%
% Input:
% param = [R Ne]
% R 	= radius (meters)
% Ne	= Number of edges = 12 x 4^n, n integer
%
% Output: see object.m
% object = topol,vertex of struct object
%
% by Juan M. Rius, December 1996, modified Jan 1998 for struct object
%

function obj = sphere_Ivan(param)

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'edges_c',[],'un',[],'ds',[],'ln',[],'cent',[]);

R = param(1); Ne = param(2);

n = log(Ne/12)/log(4);
view_var('n',n);
view_var('Mesh size',2*pi*R/(4*2^n));
if floor(n)~=n, error('Ne must be = 12 x 4^n, n integer'); end

obj.vertex = R * [ 0  1  0 -1  0  0; % Cartesian coordinates of obj.vertex
              	   0  0  1  0 -1  0;
               	   1  0  0  0  0 -1];

obj.topol = [ 1 1 1 1 6 6 6 6;
	  	      3 4 5 2 2 3 4 5;
		      2 3 4 5 3 4 5 2];

Nt = 8; Nv = 6;
for it = 1:n,			% Divide by 2 mesh size
for t = 1:Nt,			% For each triangle
	v1 = obj.topol(1,t); v2 = obj.topol(2,t); v3 = obj.topol(3,t);
	r12 = (obj.vertex(:,v1) + obj.vertex(:,v2))/2; r12 = r12*R/norm(r12);
	r23 = (obj.vertex(:,v2) + obj.vertex(:,v3))/2; r23 = r23*R/norm(r23);
	r31 = (obj.vertex(:,v3) + obj.vertex(:,v1))/2; r31 = r31*R/norm(r31);

	% Check if new obj.vertex already exists
	e12 = find(all([r12(1)==obj.vertex(1,:); r12(2)==obj.vertex(2,:); r12(3)==obj.vertex(3,:)]));
	e23 = find(all([r23(1)==obj.vertex(1,:); r23(2)==obj.vertex(2,:); r23(3)==obj.vertex(3,:)]));
	e31 = find(all([r31(1)==obj.vertex(1,:); r31(2)==obj.vertex(2,:); r31(3)==obj.vertex(3,:)]));
	if e12, v12 = e12; else v12 = Nv+1; Nv = Nv+1; obj.vertex = [obj.vertex r12]; end
	if e23, v23 = e23; else v23 = Nv+1; Nv = Nv+1; obj.vertex = [obj.vertex r23]; end
	if e31, v31 = e31; else v31 = Nv+1; Nv = Nv+1; obj.vertex = [obj.vertex r31]; end

	obj.topol(:,t) = [v12; v23; v31]; % Replace current triangle
	obj.topol = [obj.topol [v2; v23; v12] [v23; v3; v31] [v31; v1; v12]];
	end
	Nt = Nt*4;
end


if(length(param)>2)
sprintf('%s','unesi rot offset')
rot=input('rot=');  
 Rrot = [cos(rot) -sin(rot) 0; ...
      sin(rot)  cos(rot) 0; ...  % rotacijska matrica sa slit napisana
              0           0  1];
 obj.vertex = Rrot*[obj.vertex(1,:);obj.vertex(2,:);obj.vertex(3,:)];
end




