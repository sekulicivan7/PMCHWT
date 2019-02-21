% setradio(uih)
% Sets radiobutton uih to ON and the other related radiobuttons to OFF
%
% uih	= handle to checkbox or radiobutton
%
% version 3.0, Juan M. Rius, Sept 1996

function setradio(uih)

UserData = gud(uih);

for n = 2:UserData(1)+1,
	nuih = UserData(n);
	if nuih == uih, set(nuih,'Value',get(nuih,'Max'));
	else, set(nuih,'Value',get(nuih,'Min'));
	end;
end;
