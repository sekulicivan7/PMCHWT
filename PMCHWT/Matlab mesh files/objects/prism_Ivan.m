% obj = prism_Ivan(param)
% 
% Each number of basis functions on each edge of the prism
%	Nx -> param(2)
%	Ny -> param(4) 
%	Nz -> param(6)
% Each prism dimensions :
%	Lx (param(1)),
% 	Ly (param(3))
%	Lz (param(5)).
%
% by Ed. Ubeda, february 2016
%

function obj = prism_Ivan(param)

obj = struct('vertex',[],'topol',[],'trian',[],'edges',[],'edges_c',[],'un',[],'ds',[],'ln',[],'cent',[]);

Lx = param(1); Nx = param(2);
Ly = param(3); Ny = param(4);
Lz = param(5); Nz = param(6);

%primer discretitzo les arestes
arx = (linspace(-Lx/2, Lx/2, Nx+1)); %Nx funcions base
ary = (linspace(-Ly/2, Ly/2, Ny+1)); %Ny funcions base
arz = (linspace(-Lz/2, Lz/2, Nz+1)); %Nz funcions base

%Ara efectuo la discretització del prisma a partir de les arestes

%cares definides per x i y (són dues simètriques)
[c2 c1]=meshgrid(ary,arx);
zax = c1(:).';
zay = c2(:).';
zaz = -Lz/2*ones(1,(Nx+1)*(Ny+1));

 
%cares definides per y i z (són dues simètriques)
[c1 c2]=meshgrid(ary,arz);
zby = c1(:).';
zbz = c2(:).';
zbx = -Lx/2*ones(1,(Ny+1)*(Nz+1));

%cares definides per z i x (són dues simètriques)
[c1 c2]=meshgrid(arx,arz);
zcx = c1(:).';
zcz = c2(:).';
zcy = -Ly/2*ones(1,(Nz+1)*(Nx+1));

clear c1 c2 ;

% El format d'entrega de les dades són dues matrius z1 i z2:
% -----> z1 és de dimensions 3*Nt (número de triangles), essent que
%        a cada triangle se li assigna una terna que indica els números 
%        de node (vèrtex) associats a aquell triangle en qüestió.
% -----> z2 és de dimensions 3*Nv (número de vèrtexs), a on es 
%        presenten els tres components de cada vèrtex.

vec=[];
nac=Nx+3;
while nac<=(Ny-1)*(Nx+1)+2,
   vec = [vec nac:nac+Nx-2];
   nac = nac+Nx+1;
end

 
z2x = [zcx -zbx(Nz+2:(Ny+1)*(Nz+1)) zcx(Nx*(Nz+1):-1:1) zbx(Ny*(Nz+1):-1:Nz+2) zax(vec) zax(vec)];
z2y = [zcy zby(Nz+2:(Ny+1)*(Nz+1)) -zcy(Nx*(Nz+1):-1:1) zby(Ny*(Nz+1):-1:Nz+2) zay(vec) zay(vec)];
z2z = [zcz zbz(Nz+2:(Ny+1)*(Nz+1)) zcz(Nx*(Nz+1):-1:1) zbz(Ny*(Nz+1):-1:Nz+2) -zaz(vec) zaz(vec)];

if(length(param)>6)
sprintf('%s','unesi x offset')
xpom=input('xpom=');  
sprintf('%s','unesi y offset')
ypom=input('ypom=');
sprintf('%s','unesi z offset')
zpom=input('zpom=');
obj.vertex = [z2x+xpom;z2y+ypom;z2z+zpom];
else
obj.vertex = [z2x;z2y;z2z];    
end
clear z2x z2y z2z;

% A continuació obtinc la matriu de triangles z1 (ztop)
% -->Triangles de les cares zx i zy primeres de la rotació
z11=[];
for n = 1:Nx+Ny,
 for r = 1:Nz,
  aux = [(n-1)*(Nz+1)+r n*(Nz+1)+r+1 n*(Nz+1)+r;
	 (n-1)*(Nz+1)+r (n-1)*(Nz+1)+r+1 n*(Nz+1)+r+1];
  z11 = [z11; aux];
 end
end

% Numero la filera on hi ha el canvi de sentit
z12=[];
for r = 1:Nz,
	if r==1,
		aux=[(Nx+Ny+1)*(Nz+1)-r (Nx+Ny+1)*(Nz+1)+r-1 (Nx+Ny+1)*(Nz+1)+r;
			(Nx+Ny+1)*(Nz+1)-r (Nx+Ny+1)*(Nz+1)+r (Nx+Ny+1)*(Nz+1)+r+1];
	else
		aux=[(Nx+Ny+1)*(Nz+1)-r (Nx+Ny+1)*(Nz+1)-r+1 (Nx+Ny+1)*(Nz+1)+r;
			(Nx+Ny+1)*(Nz+1)-r (Nx+Ny+1)*(Nz+1)+r (Nx+Ny+1)*(Nz+1)+r+1];
	end
	z12=[z12; aux];
end

% M'encarrego de la resta de parets laterals, ara ja amb el sentit de
% numeració ja canviat
z13=[];
	for n=1:Nx+Ny-2,
		for r=1:Nz,
			aux=[(Nx+Ny+1)*(Nz+1)+(n-1)*(Nz+1)+r (Nx+Ny+1)*(Nz+1)+n*(Nz+1)+r (Nx+Ny+1)*(Nz+1)+(n-1)*(Nz+1)+r+1;
				(Nx+Ny+1)*(Nz+1)+(n-1)*(Nz+1)+r+1 (Nx+Ny+1)*(Nz+1)+n*(Nz+1)+r  (Nx+Ny+1)*(Nz+1)+n*(Nz+1)+r+1];  
			z13=[z13; aux];
		end 
	end   

%Per últim, la darrera filera, que tanca les parets laterals
z14=[];
	for r=1:Nz,
		aux=[(2*(Nx+Ny)-1)*(Nz+1)+r Nz+2-r (2*(Nx+Ny)-1)*(Nz+1)+r+1;
			(2*(Nx+Ny)-1)*(Nz+1)+r+1 Nz+2-r Nz+2-r-1];
		z14=[z14; aux];
	end



%tapa superior (z=+Lz/2)
z15=[]; 
%primera filera

	for r=1:Nx,
		if (r==1),
			aux = [Nz+1 (2*(Nx+Ny)-1)*(Nz+1)+1 2*(Nz+1);
				(2*(Nx+Ny)-1)*(Nz+1)+1 2*(Nx+Ny)*(Nz+1)+1 2*(Nz+1)];
		end

		if (r>1)&(r<Nx),
			aux = [r*(Nz+1) 2*(Nx+Ny)*(Nz+1)+r-1 (r+1)*(Nz+1);
				2*(Nx+Ny)*(Nz+1)+r-1 2*(Nx+Ny)*(Nz+1)+r (r+1)*(Nz+1)];
		end

		if (r==Nx),
			aux = [Nx*(Nz+1) 2*(Nx+Ny)*(Nz+1)+Nx-1 (Nx+1)*(Nz+1);
				2*(Nx+Ny)*(Nz+1)+Nx-1 (Nx+2)*(Nz+1) (Nx+1)*(Nz+1)];
		end
		z15 = [z15; aux];
	end

%segona filera
%part de dalt
    for n=1:Ny-2,
     aux = [(2*(Nx+Ny)-n)*(Nz+1)+1 (2*(Nx+Ny)-(n+1))*(Nz+1)+1 ...
             2*(Nx+Ny)*(Nz+1)+(n-1)*(Nx-1)+1;
	    (2*(Nx+Ny)-(n+1))*(Nz+1)+1 ...
	     2*(Nx+Ny)*(Nz+1)+n*(Nx-1)+1 ...
            2*(Nx+Ny)*(Nz+1)+(n-1)*(Nx-1)+1];
      z15 = [z15; aux];
    end
  
    
   %part d'enmig
      for n=1:Nx-2,
       for r=1:Ny-2,
        aux = [2*(Nx+Ny)*(Nz+1)+(r-1)*(Nx-1)+n ...
	       2*(Nx+Ny)*(Nz+1)+r*(Nx-1)+n ...
	       2*(Nx+Ny)*(Nz+1)+(r-1)*(Nx-1)+n+1;
	       2*(Nx+Ny)*(Nz+1)+r*(Nx-1)+n ...
	       2*(Nx+Ny)*(Nz+1)+r*(Nx-1)+n+1 ...
	       2*(Nx+Ny)*(Nz+1)+(r-1)*(Nx-1)+n+1];
        z15 = [z15; aux];
       end
      end
    

   %part d'abaix
      for n=1:Ny-2,
        aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*n ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(n+1) ...
	       (Nx+1+n)*(Nz+1);
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(n+1) ...
	       (Nx+1+n+1)*(Nz+1) ...
	       (Nx+1+n)*(Nz+1)];
        z15 = [z15; aux];
      end	       
   

 %tercera filera
 for r=1:Nx,
  if r==1,
   aux = [(Nx+Ny)*(Nz+1) 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1) ...
          (Nx+Ny+1)*(Nz+1)+1;
	  (Nx+Ny)*(Nz+1) (Nx+Ny+1)*(Nz+1)+1 (Nx+Ny+1)*(Nz+1)];
  end

  if (r>1)&(r<Nx),
   aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)-(r-2) ...
          2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)-(r-1) (Nx+Ny+1+(r-1))*(Nz+1)+1;
	  2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)-(r-2) ...
          (Nx+Ny+1+(r-1))*(Nz+1)+1 (Nx+Ny+1+(r-2))*(Nz+1)+1];
  end

  if r==Nx,
   aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-2)+1 ...
	  (2*Nx+Ny+1)*(Nz+1)+1 (2*Nx+Ny)*(Nz+1)+1;
	  2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-2)+1 ...
	  (2*Nx+Ny)*(Nz+1)+1 (2*Nx+Ny-1)*(Nz+1)+1];
  end
  z15 = [z15; aux];

 end


%tapa inferior (z=-Lz/2)
z16=[]; 
 %primera filera

 for r=1:Nx,
  if r==1,
   aux = [(Nz+1)+1 2*(Nx+Ny)*(Nz+1) 1;
	  (Nz+1)+1 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+1 2*(Nx+Ny)*(Nz+1)];
  end

  if (r>1)&(r<Nx),
   aux = [r*(Nz+1)+1 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r-1 (r-1)*(Nz+1)+1 ;	
	  r*(Nz+1)+1 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r ...
		2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r-1 ];
  end

  if r==Nx,
   aux = [Nx*(Nz+1)+1 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+Nx-1 ...
	  (Nx-1)*(Nz+1)+1;
	  Nx*(Nz+1)+1 (Nx+1)*(Nz+1)+1 2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+Nx-1];
  end
  z16 = [z16; aux];
 end

 %segona filera
   %part de dalt
    for n=1:Ny-2,
     aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(n-1)*(Nx-1)+1 ...
	    (2*(Nx+Ny)-n)*(Nz+1) (2*(Nx+Ny)-(n-1))*(Nz+1);
	    2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(n-1)*(Nx-1)+1 ...
	    2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+n*(Nx-1)+1 ...
	    (2*(Nx+Ny)-n)*(Nz+1)];
      z16 = [z16; aux];
    end
    
   %part d'enmig
      for n=1:Nx-2,
       for r=1:Ny-2,
        aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(r-1)*(Nx-1)+n ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(r-1)*(Nx-1)+n+1 ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r*(Nx-1)+n;
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r*(Nx-1)+n ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(r-1)*(Nx-1)+n+1 ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+r*(Nx-1)+n+1];
        z16 = [z16; aux];
       end
      end

   %part d'abaix
      for n=1:Ny-2,
        aux = [(Nx+1+n-1)*(Nz+1)+1 ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(n+1)*(Nx-1) ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+n*(Nx-1) ;
	       (Nx+1+n-1)*(Nz+1)+1  ...
	       (Nx+1+n)*(Nz+1)+1  ...
	       2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(n+1)*(Nx-1)];
        z16 = [z16; aux];
      end	


%tercera filera
for r=1:Nx,
	if r==1,
		aux = [(Nx+Ny-1)*(Nz+1)+1 (Nx+Ny+2)*(Nz+1) ...
			2*(Nx+Ny)*(Nz+1)+2*(Nx-1)*(Ny-1);
			(Nx+Ny-1)*(Nz+1)+1 (Nx+Ny)*(Nz+1)+1 ...
			(Nx+Ny+2)*(Nz+1)];
	end

	if (r>1)&(r<Nx),
		aux = [2*(Nx+Ny)*(Nz+1)+2*(Nx-1)*(Ny-1)-(r-2) ...
			(Nx+Ny+2+(r-1))*(Nz+1) 2*(Nx+Ny)*(Nz+1)+2*(Nx-1)*(Ny-1)-(r-1);
			2*(Nx+Ny)*(Nz+1)+2*(Nx-1)*(Ny-1)-(r-2) ...
			(Nx+Ny+2+(r-2))*(Nz+1) (Nx+Ny+2+(r-1))*(Nz+1)];
	end

	if r==Nx,
		aux = [2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(Nx-1)*(Ny-2)+1 ...
			(2*Nx+Ny+1)*(Nz+1) (2*Nx+Ny+2)*(Nz+1);
			2*(Nx+Ny)*(Nz+1)+(Nx-1)*(Ny-1)+(Nx-1)*(Ny-2)+1 ...
		  	(2*Nx+Ny)*(Nz+1) (2*Nx+Ny+1)*(Nz+1)];
	end
	z16 = [z16; aux];

end   

zlat =  [z11; z12; z13; z14];
Nlat = length( zlat );

obj.topol = [zlat; z15; z16].';
