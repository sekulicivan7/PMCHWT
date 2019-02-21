% [namelist,var1mat,var2mat] = savename(varname,name,namelist,var1,var1mat,var2,var2mat)
% Save results in list. If new name matches with an already saved name,
% a dialog window asks: replace, append (at end), skip
%
% Input:
% varname  = name of variable (string)
% name	   = name of new entry to save
% namelist = list of saved names
% var1	   = variable to save (vector)
% var1mat  = matrix of saved variables (vectros in columns)
% var2	   = 2nd variable to save (optional)
% var2mat  = 2nd matrix of saved variables (optional)
%
% Output:
% Updated variables
%
% ver 3.3, Juan M. Rius, Jan 1997


function [namelist,var1mat,var2mat] = savename(varname,name,namelist,var1,var1mat,var2,var2mat)

if ~length(var1), return;	% Nothing to save
end

inlist = findname(namelist,name);

if all(~inlist),	% If name does no exist in list
	if length(namelist), 	namelist = str2mat(namelist,name);
	else,			namelist = name;
	end
	var1mat = [var1mat var1];
	if nargin > 5, var2mat = [var2mat var2];
	end
else,	% for necessary if name matches with more than one entry
	for i = find(inlist)',
		text = [varname ': ''' name ''' matches with saved result ''' deblank(namelist(i,:)) ''''];
		clk = questdlg(text,'Atention:','Replace','Append','Skip','Skip');
      
      switch clk,
      case 'Replace',
         [m,n] = size(namelist);
         if m==1, namelist = name;
         else
            switch i,
            case m, namelist = str2mat(namelist(1:i-1,:),name);
            case 1, namelist = str2mat(name,namelist(i+1:m,:));
            otherwise, namelist = str2mat(namelist(1:i-1,:),name,namelist(i+1:m,:));
            end
         end
         var1mat(:,i) = var1;
         if nargin > 5, var2mat(:,i) = var2;
         end
		case 'Append',
			if length(namelist), 	namelist = str2mat(namelist,name);
			else,			namelist = name;
			end
			var1mat = [var1mat var1];
			if nargin > 5, var2mat = [var2mat var2];
			end
		end
	end
end;
