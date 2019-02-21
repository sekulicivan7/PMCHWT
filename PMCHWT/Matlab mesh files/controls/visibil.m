% visibil(uih)
% Calls updatvis to update visibility of handles that depend on uih
%
% uih	= handle to checkbox or radiobutton array
%	  if uih is a vector, only the 1st element is considered,
%	  so it must be a checkbutton array
%
% version 3.0, Juan M. Rius, Sept 1996

function visibil(uih)

if length(uih)>1, uih = uih(1);
end

mat = gud(uih);

if	strcmp(get(uih,'Style'),'checkbox'),	n=1;
elseif	strcmp(get(uih,'Style'),'radiobutton'),	n=mat(1)+2;
else 	error('Handle must be radiobutton array or checkbox');
end

updatvis(mat(n:length(mat)));

