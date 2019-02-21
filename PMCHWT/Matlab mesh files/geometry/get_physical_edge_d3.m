% Procedure to obtain the indices for the physical and non-physical edges
%
% Physical edges = edges shaped by two triangles that are non-coplanar
%
% Output parameter:
% set_phys = set of physical edges
% set_non_phys = set of non-physical edges
%
% by Ed. Ubeda, September 2012

function [ set_phys, set_non_phys,set_physT] = get_physical_edge_d3(obj);

Ne = length( obj.ln );

set_phys = [];
set_non_phys = [];
set_physTp=[];
set_physTm=[];

for n=1:Ne,
    
    Tp = obj.edges(1,n);
    Tm = obj.edges(2,n);
    
    diff = abs( norm( obj.un(:,Tp) - obj.un(:,Tm) ) );
    
    if (diff<1000*eps),
        set_non_phys = [set_non_phys n ];
    elseif (diff>=1000*eps),
        set_phys = [set_phys n ];
        set_physTp=[set_physTp,Tp];
        set_physTm=[set_physTm,Tm];
    end;  %%% if (diff<10*eps),
    
end; %%% for p=1:Ne,

set_physT=unique([set_physTp,set_physTm]);
