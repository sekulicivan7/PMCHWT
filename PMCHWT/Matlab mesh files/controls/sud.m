% sud(handle,value,n,m)
% Store User Data in graphics object
% Useful for GUI's that store variables in UserData arrays
%
% handle = handle to graphics object
% value  = value to store
% n	 = row index to UserData matrix (optional)
% m	 = column index to UserData matrix (optional)
%	 n or m ==0 are equivalent to : (whole row or column)
%
% version 3.1, Juan M. Rius, Oct 1996

function sud(handle,value,n,m)

if nargin==2,		set(handle,'UserData',value);
elseif nargin==3,	UserData = get(handle,'UserData');
			N = length(UserData);
			if ~n, n = 1:N; end;
			UserData(n) = value;
			set(handle,'UserData',UserData);
elseif nargin==4,	UserData = get(handle,'UserData');
			[N,M] = size(UserData);
			if ~n, n = 1:N; end;
			if ~m, m = 1:M; end;
			UserData(n,m) = value;
			set(handle,'UserData',UserData);
else error('Incorrect argument count');
end



