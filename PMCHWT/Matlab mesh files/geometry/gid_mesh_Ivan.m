% obj = gid_mesh(p_gid)
%
% Input:
% p_gid: struct with the following field
%		file = GID mesh file
%
% Output: see object.m
% obj = object struct, with fields vertex,topol and Ng
%
% Notes:
% Node numbers must begin with 1 and be correlative
% Mesh must contain only triangular elements
% The vertex matrix begins after the keyword 'Coordinates'
% The topology matrix begins after the keyword 'Elements'
% 
% Juan M. Rius, Josep Parron June 1999
%

function obj = gid_mesh_Ivan(p_gid)

fich = p_gid.file;
if(iscell(fich)) fich = fich{1};
end

% For compatibility with present version of run_3d
if(~ischar(fich)) error('Argument fich must be string or string cell');
end

fid = fopen(fich);
if(fid==-1) error(sprintf('Cannot open file %s',fich));
end

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'un',[],'ds',[],'ln',[],'cent',[]);

tmp = fgetl(fid);
while ~strcmp(tmp,'Coordinates');
   tmp = fgetl(fid);
end
tmp = fscanf(fid,'%f',inf);
tmp = reshape(tmp,4,length(tmp)/4);
if( any( tmp(1,:)~=(1:size(tmp,2)) ) ) error('Node numbers are not correlative');
end
obj.vertex = tmp(2:4,:);


tmp = fgetl(fid);
while ~strcmp(tmp,'Elements');
   tmp = fgetl(fid);
end
tmp = fscanf(fid,'%f',inf);
tmp = reshape(tmp,4,length(tmp)/4);
obj.topol = tmp(2:4,:);

fclose(fid);

obj = get_edge_Ivan(obj);


   

