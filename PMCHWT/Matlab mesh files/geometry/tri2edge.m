
function set_edge=tri2edge(obj)

% reordering of the matrix to edge-wise
Ne = size(obj.edges,2);

for n=1:Ne
    
    TP=obj.edges(1,n);
    TM=obj.edges(2,n);
    
    i_p= find(abs(obj.trian(:,TP))==n);
    i_m= find(abs(obj.trian(:,TM))==n);
    
    set_edge(n) = 3*(TP-1) + i_p;
    set_edge(n+Ne) = 3*(TM-1) + i_m;
    
end