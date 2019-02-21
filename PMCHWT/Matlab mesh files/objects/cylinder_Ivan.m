% A cylinder is returned
%
% param(1)  -----> cylinder radius 
% param(2)  -----> steps on the R direction 
% param(3)  -----> cylinder trunk height 
% param(4)  -----> steps on the z direction 
% param(5)  -----> steps on the phi direction (at the highest level) 
%
%
% by Eduard Ubeda Farre, December 1998.
%

function obj = cylinder_Ivan(param)

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'edges_c',[],'un',[],'ds',[],'ln',[],'cent',[]);

Rd = param(1);   % Cylinder Radius
Nrd = param(2);  

Lz = param(3); 	% Cylinder height
Nz = param(4);  

Npol = param(5);  % Steps on the Phi direction

arR = linspace(0, Rd, Nrd+1); 	       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the TOP cylinder BASIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj.vertex = [0 0 0].';

set0 = 1;
Nphi0 = 1;

for n=2:(Nrd+1),
   
   Nph_int = (n-1)*Npol;			% number of vertex inside of the row
   
   % obj.vert update
   
   auxphi = linspace(0, 2*pi , Nph_int+1);
   auxphi = auxphi(1:Nph_int);
   auxR = arR(n)*ones(1,Nph_int);
   auxZ = zeros(1,Nph_int);
   
   % Transformation from cylindrical coordinates into cartesian coordinates
   
   auxX = auxR .*  cos(auxphi);
   auxY = auxR .*  sin(auxphi);
   
   obj.vertex = [ obj.vertex [auxX; auxY; auxZ] ];
   
   % obj.topol update
   
   Nphiff = size( obj.vertex , 2);
   
   setf = (set0(length(set0))+1):Nphiff;
   Nphif = length(setf);
   
   auxT = [];
   
   for m=1:Npol,
      
      for q=1:(n-1),
         
         auxT = [auxT  ...
               [ setf( rem( (m-1)* (n-1) + q - 1 , Nphif ) + 1 ) ...
                  set0( rem( (m-1)*(n-2) + q - 1 , Nphi0) + 1 ) ...
                  setf( rem( (m-1)* (n-1) + q + 1 - 1 , Nphif) + 1 ) ...
               ].' ];
         
         if q~=(n-1),		
            
            auxT = [auxT ...
                  [ set0( rem( (m-1)*(n-2) + q - 1 , Nphi0) + 1 ) ...
                     set0( rem( (m-1)*(n-2) + q + 1 - 1 , Nphi0 ) + 1 )  ...
                     setf( rem( (m-1)*(n-1) + q + 1 - 1, Nphif ) + 1 )  ...
                  ].' ];            
            
         end  %% if q~=(n-1),
         
      end  %% for q=1:n-1,

   end  %% for m=1:Npol
   
   obj.topol = [obj.topol auxT];
   
   set0 = setf;
   Nphi0 = Nphif;
   
end	%% for n=2:Nrd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the cilinder TRUNK %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% obj.vertex update

Npol2 = Nrd*Npol;

arphi=linspace(0,2*pi,Npol2+1).';
z=linspace(0,-Lz,Nz+1);
z=z(2:length(z));

arX = Rd * cos(arphi(1:Npol2)) * ones( size(z) );
arY = Rd * sin(arphi(1:Npol2)) * ones( size(z) );
arZ = ones(Npol2,1)*z;

X = arX(:).';
Y = arY(:).';
Z = arZ(:).';

N0 = size(obj.vertex,2);

obj.vertex=[obj.vertex [X; Y; Z]];

%% obj.topol update

for n=1:Nz,
   for r=1:Npol2,
      if r~=Npol2,
         aux=[N0+(n-2)*Npol2+r   N0+(n-2)*Npol2+r+1    N0+(n-1)*Npol2+r      ;
            N0+(n-2)*Npol2+r+1   N0+(n-1)*Npol2+r+1    N0+(n-1)*Npol2+r      ]; 
      else,
         aux=[N0+(n-2)*Npol2+Npol2  N0+(n-2)*Npol2+1   N0+(n-1)*Npol2+Npol2   ;
            N0+(n-2)*Npol2+1        N0+(n-1)*Npol2+1   N0+(n-1)*Npol2+Npol2   ]; 
      end
      obj.topol=[obj.topol aux.'];
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the BOTTOM cylinder BASIS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nor = size(obj.vertex,2);
set0 = (Nor-Npol2+1):Nor;			%% most exterior vertex row
Nphi0 = length( set0 );

for n=Nrd:-1:2,
   
   Nph_int = (n-1)*Npol;			% number of vertex inside of the row
   
   % obj.vert update
   
   auxphi = linspace(0, 2*pi, Nph_int+1);
   auxphi = auxphi(1:Nph_int);
   auxR = arR(n)*ones(1,Nph_int);
   auxZ = -Lz*ones(1,Nph_int);
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
   for m=1:Npol,
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
      end  %% for q=1:(n-1),
      
   end  %% for m=1:Npol
   %
   obj.topol = [obj.topol auxT];
   
   set0 = setf;
   Nphi0 = Nphif;
   %
end	%% for n=1:Nr

%% final update 

obj.vertex = [obj.vertex [0 0 -Lz].'];
Nv = size(obj.vertex,2);

auxT=[];
for m=1:Npol,
   auxT = [auxT [set0( rem(m-1,Nphi0)+1 ) set0( rem(m+1-1,Nphi0)+1 )  Nv ].' ];
end

obj.topol = [obj.topol auxT];


