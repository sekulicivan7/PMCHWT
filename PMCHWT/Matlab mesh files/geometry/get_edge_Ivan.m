% object_out = get_edge(object_in)
% Computes edges, ln, trian, un, ds and cent fields of object
%
%   field==1       --> EFIE[R]  [ed]
%   field==2       --> MFIE[R]  [ed]
%   field==3       --> MFIE[monoR]  [ed]
%   field==4       --> EFIE[tet][monoR]  [tri]
%   field==5       --> EFIE[wed][monoR]  [tri]
%   field==6       --> EFIE[tet][eoR]  [ed]
%   field==7       --> EFIE[wed][eoR]  [ed]
%   field==8       --> EFIE[TN][monoR] [tri]
%
% by Ed. Ubeda [February 2016]

function obj = get_edge_Ivan(object_in)

obj = struct('vertex',object_in.vertex,'topol',object_in.topol,'trian',[],'edges',[],'un',[],'ds',[],'ln',[],'cent',[]);

v1 = obj.vertex(:,obj.topol(1,:));
v2 = obj.vertex(:,obj.topol(2,:));
v3 = obj.vertex(:,obj.topol(3,:));

obj.cent =(v1+v2+v3)/3;
c    = cross(v3-v1,v2-v1);
obj.un   = unitary(c);
obj.ds   = sqrt(sum(c.^2))/2;

Nt = length(obj.ds);
obj.trian = zeros(3,Nt);
obj.edges = zeros(4,ceil(Nt*3/2));

eg = 1;								% Global number of current edge
for Tp = 1:Nt,						% Triangle T+
    for el = [3 2 1],					% Local edge, scan order [3 2 1] for compatibility
        if ~obj.trian(el,Tp),		% Edge not found yet
            % Find vertices of this edge
            ver = obj.topol(find([1 2 3]~=el),Tp);
            if el == 2, ver = flipud(ver);
            end;
            
            [tmp T1] = find(obj.topol == ver(1)); % T1 = triangles that have ver(1)
            [tmp T2] = find(obj.topol == ver(2)); % T2 = triangles that have ver(2)
            
            [TT1 TT2] = meshgrid(T1,T2);
            Tcom = TT1(TT1==TT2);	% Triangles with common edge
            
            if length(Tcom) == 2,
                obj.trian(el,Tp) = eg;
                obj.edges(1,eg) = Tp;		% T+
                obj.edges(3,eg) = obj.topol(el,Tp);
                
                Tm = Tcom(find(Tcom~=Tp)); 	% T-, not equal to Tp
                obj.edges(2,eg) = Tm;
                
                % Find the local number of eg in T-
                v = find(obj.topol(:,Tm)~=ver(1) & (obj.topol(:,Tm)~=ver(2)));
                obj.trian(v,Tm) = -eg;
                obj.edges(4,eg) = obj.topol(v,Tm);
                
                obj.ln(eg) = norm(obj.vertex(:,ver(1))-obj.vertex(:,ver(2)));
                eg = eg+1;
            elseif length(Tcom) > 2, 
                error('More that 2 triangles share an edge');
            end;
        end;
    end;
end;

obj.edges = obj.edges(:,1:eg-1);   % Remove void edges, in case of open object






