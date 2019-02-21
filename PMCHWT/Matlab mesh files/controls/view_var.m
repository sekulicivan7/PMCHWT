% view_var(text,x)
% View nicely the contents of variable with prefix text
%
% Version 1.0, Oct. 1996, Juan M. Rius


function view_var(text,x)

st = [text, ' = '];

if isstr(x),
	disp([st x]);
	return
end

[N,M] = size(x);

stl = length(st);

if N*M>1,
	st = [st, '['];
	maxel = zeros(1,M-1);
	for n=1:N,
		for m=1:M-1,
			tmp = length(num2str(x(n,m)));
			if tmp > maxel(m), maxel(m) = tmp; end
		end
	end

	for n=1:N,
		for m=1:M-1,
			stel = num2str(x(n,m));
			st = [st, ' ', stel];

			for i=length(stel):maxel(m),
				st = [st, ' '];
			end
		end

		st = [st, ' ', num2str(x(n,M))];
		st = [st, sprintf('\n')];

		if n~=N, for i=1:stl+1,
			st = [st, ' '];
		end; end
	end
	st = [st(1:length(st)-1), ' ]'];
else
	st = [st, num2str(x)];
end
disp(st);

% num2str is necessary for displaying complex numbers
