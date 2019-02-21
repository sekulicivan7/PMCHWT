% UserData =  gud(handle,n,m)
% Get User Data from graphics object
% Useful for GUI's that store variables in UserData arrays
%
% handle = handle to graphics object
% n	 = row index to UserData matrix (optional)
% m	 = column index to UserData matrix (optional)
%	 n or m ==0 are equivalent to : (whole row or column)
%
% Returns UserData or UserData(index)
%
% version 3.1, Juan M. Rius, Oct 1996

function UserData =  gud(handle,n,m)

UserData =get(handle,'UserData');

% if nargin==2,		N = length(UserData);
% 			if ~n, n = 1:N;
%             n
%             end;
% 			UserData = UserData(n)
%             
% elseif nargin==3,	[N,M] = size(UserData);
% 			if ~n, n = 1:N; 
%             n
%             end;
% 			if ~m, m = 1:M; 
%             m
%             end;
% 			UserData = UserData(n,m)
            
end


