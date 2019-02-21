% A cone is returned
%
% param(1)  -----> cone height
% param(2)  -----> cone radius
% param(3)  -----> steps on the z direction
% param(4)  -----> steps on the phi direction (at the highest level)
%
%  
% by Eduard Ubeda Farre, February 2016.
%


function obj = cone_Ivan(param)

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'edges_c',[],'un',[],'ds',[],'ln',[],'cent',[]);

Lz = param(1);
Rd = param(2);

Nz = param(3);
Nphi00 = param(4);

arz = linspace(Lz, 0 ,Nz+1); 	       %Nz funcions base
arR = linspace(0, Rd, Nz+1); 	       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the CONE STRUCTURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj.vertex = [ 0 0 Lz].';

for n=1:Nz,
	%
	Nphi0 = size(obj.vertex,2);	%vertex at the origin
	Nphi = n*Nphi00;		%number of vertex of this row (n vertex per sector)
	Nphi1 = (n-1)*Nphi00;		%number of vertex of the precedent
					% row ((n-1) vertex per sector)
	%
	% obj.vert update
	%
	auxphi = linspace( 0 , 2*pi, Nphi+1);   %Nphi funcions base
	auxphi = auxphi(1:Nphi);
	auxR = arR(n+1)*ones(1,Nphi);
	auxZ = arz(n+1)*ones(1,Nphi);
	%
	% transformation from cilindrical coordinates into cartesian coordinates
	%
	auxX = auxR .*  cos(auxphi);
	auxY = auxR .*  sin(auxphi);
	%
	obj.vertex = [ obj.vertex [auxX; auxY; auxZ] ];
	%
	% obj.topol update
	%
	Nphif = size( obj.vertex , 2);
	%
	set0 = (Nphi0-Nphi1+1):Nphi0;
	setf = (Nphi0+1):Nphif;
	%
	auxT = [];
	%
	if n==1,
		for m=1:Nphi00,
			auxT = [auxT [rem(m-1,Nphi)+2 1 rem(m+1 -1,Nphi)+2].' ];
		end
    else
		for m=1:Nphi00,
			%
			for q=1:n
				auxT = [auxT  ...
				[ setf( rem( (m-1)*n + q - 1 , Nphi  ) + 1 ) ...
				  set0( rem( (m-1)*(n-1) + q - 1 , Nphi1) + 1 ) ...
				  setf( rem( (m-1)*n + q + 1 - 1 , Nphi) + 1 ) ].' ];
				%
				if q~=n,		
					%
					auxT = [auxT ...
					[ set0( rem( (m-1)*(n-1) + q - 1 , Nphi1) + 1 )  ...
					  set0( rem( (m-1)*(n-1) + q + 1 - 1 , Nphi1 ) + 1 ) ...
					  setf( rem( (m-1)*n + q + 1 - 1, Nphi ) + 1 ) ].' ];
					%
				end  %% if q~=n,
				%
			end  %% for q=1:n
			%
		end  %% for m=1:Nphi00
		%	
	end  %%else if n==1
	%
	obj.topol = [obj.topol auxT];
	%
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Build the CONE BASIS %
%%%%%%%%%%%%%%%%%%%%%%%%

set0 = setf;		%% most exterior vertex row
Nphi0 = length( set0 );
%
for n=Nz:-1:2,
	%
	Nph_int = (n-1)*Nphi00;			% number of vertex inside of the row
	%
	% obj.vert update
	%
	auxphi = linspace(0, 2*pi , Nph_int+1);
	auxphi = auxphi(1:Nph_int);
	auxR = arR(n)*ones(1,Nph_int);
	auxZ = zeros(1,Nph_int);
	%
	% transformation from cilindrical coordinates into cartesian coordinates
	%
	auxX = auxR .*  cos(auxphi);
	auxY = auxR .*  sin(auxphi);
	%
	obj.vertex = [ obj.vertex [auxX; auxY; auxZ] ];
	%
	% obj.topol update
	%
	Nphiff = size( obj.vertex , 2);
	%
	setf = (set0(length(set0))+1):Nphiff;
	Nphif = length(setf);
	%
	auxT = [];
	%
	for m=1:Nphi00,
		%
		for q=1:n,
			auxT = [auxT  ...
			[ set0( rem( (m-1)*n + q - 1 , Nphi0  ) + 1 ) ...
			  set0( rem( (m-1)*n + q + 1 - 1 , Nphi0) + 1 ) ...
			  setf( rem( (m-1)*(n-1) + q - 1 , Nphif) + 1 ) ].' ];
			%
			if q~=n,		
				%
				auxT = [auxT ...
				[ setf( rem( (m-1)*(n-1) + q - 1 , Nphif) + 1 ) ...
				  set0( rem( (m-1)*n + q + 1 - 1, Nphi0 ) + 1 ) ...
				  setf( rem( (m-1)*(n-1) + q + 1 - 1 , Nphif ) + 1 ) ].' ];
				%
			end  %% if q~=n,
			%
		end  %% for q=1:n,
		%
	end  %% for m=1:Nphi00
	%	
	obj.topol = [obj.topol auxT];
	%
	set0 = setf;
	Nphi0 = Nphif;
	%
end	%% for n=1:Nr

%% final update 

obj.vertex = [obj.vertex [0 0 0].'];

Nv = size(obj.vertex,2);


auxT=[];
for m=1:Nphi00,
	auxT = [auxT [set0( rem(m-1,Nphi0)+1 ) set0( rem(m+1-1,Nphi0)+1 )  Nv ].' ];
end

obj.topol = [obj.topol auxT];