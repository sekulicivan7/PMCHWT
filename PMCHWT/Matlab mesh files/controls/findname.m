% indb = findname(strsave,name)
% Find rows in string matrix that match with name
%
% Input:
% namelist = Matrix of row strings
% name	   = Name to match in namelist
%
% Output:
% indb	   = Vector of 1/0 indicating the rows of namelist
%	     that begin with name
%
% ver 3.3, Juan M. Rius, Jan 1997

function indb = findname(namelist,name)

N = size(namelist,1);

if N,	indb = zeros(N,1);
	for i=1:N,
		indb(i) = any(findstr(namelist(i,:),name));
	end
else,	indb = [];	% namelist is empty
end	

