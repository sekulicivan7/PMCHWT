% updatvis(handles)
%
% Updates visibility of handles according to current value
% of all radiobuttons and checkboxes in current figure
%
% handles =	Vector of handles to graphic objects.
%		They are processed independently.
%		Skips flags (integer numbers) between handles
%
% version 3.0, Juan M. Rius, Sept 1996

function updatvis(handles)

uicon =  [findobj(gcf,'Style','radiobutton'); findobj(gcf,'Style','checkbox')];

for i=1:length(handles),

handle = handles(i);
if rem(handle,1),	% If it is not a flag
vis = 1;
radiovec = [];		% Radiobuttons that have already been checked

for j=1:length(uicon),
uih = uicon(j);

% If this control has not been checked before
if length(radiovec)==0, enterloop = 1;
else enterloop = ~any(find(uih==radiovec));
end

if enterloop,
mat = gud(uih);
if length(mat)>0 & length(handle)>0, % Necessary for mat==handle
if any(find(mat==handle)),	% If handle depends on this control
	if strcmp(get(uih,'Style'),'checkbox'),
		n=1;
		value = get(uih,'Value');
	elseif strcmp(get(uih,'Style'),'radiobutton'),
		n=mat(1)+2;
		radiovec = [radiovec mat(2:n-1)]; % Add this radiobuttons to already-checked
		for i=1:mat(1),			  % Find the value of this radiobutton array
			uii = mat(i+1);
			value = get(uii,'Value');
			if value==get(uii,'Max'), break; % This button is ON
			end
		end
	else 	error('Handle must be radiobutton array or checkbox');
	end

	while n < length(mat),
		m1 = n+1;  while  rem(mat(m1),1), m1 = m1+1; end	% m1 -> 1st flag
		m2 = m1+1;
		if m2<=length(mat),
			while ~rem(mat(m2),1), m2 = m2+1; if m2>length(mat), break; end
			end	% m2 -> 1st handle
		end
		if any(mat(n:m1-1)==handle),			% handle is in this group
			vis = vis & any(mat(m1:m2-1)==value);	% Visible
			break;					% Handle has been found -> exit while
		end
		n = m2;
	end

end;	% if handle is in mat
end;  % if handle and mat are not empty
end;	% if uih is not in radiovec (enterloop)
end;	% for uih

if vis,	set(handle,'Visible','on');
else,	set(handle,'Visible','off');
end

end;	% if handle is not flag
end;	% for handle
