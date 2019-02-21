% un = unitari(x)
% Unit vector in the direction of vector x
% Arrays of vectors, 3 x N

function un = unitary(x)
a_norm =sqrt(sum(x.^2));
un =[x(1,:)./a_norm;
     x(2,:)./a_norm;
     x(3,:)./a_norm];
