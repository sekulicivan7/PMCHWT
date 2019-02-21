% script to read gid_mesh
%
% detects automatically for RWG or wire grid
%
% Diel_pec3 update 

fich_gidmesh = p_obj.gid.file; if(iscell(fich_gidmesh)) fich_gidmesh = fich_gidmesh{1}; end
if(~ischar(fich_gidmesh)) error('Argument fich must be string or string cell'); end
fid_gidmesh = fopen(fich_gidmesh);
if(fid_gidmesh==-1) error(sprintf('Cannot open file %s',fich_gidmesh)); end

tmp_gidmesh = fgetl(fid_gidmesh);
patch_gidmesh = findstr(tmp_gidmesh,'Triangle');
wire_gidmesh = findstr(tmp_gidmesh,'Linear');
fclose(fid_gidmesh);

if patch_gidmesh & wire_gidmesh, error('file should contain only wires or only patches'), end
if not(patch_gidmesh | wire_gidmesh), error('file should contain wires or patches'); end
if patch_gidmesh, obj = gid_mesh_Ivan(p_obj.gid); end
if wire_gidmesh, wobj=gid_wire(p_obj.gid); end

clear fich_gidmesh tmp_gidmesh fid_gidmesh patch_gidmesh wire_gidmesh 


