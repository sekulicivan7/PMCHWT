% [var1,var2] = loadname(varname,name,namelist,var1mat,var2mat)
% Load results from list. If name matches a name in list,
% a dialog window prompts for confirmation
%
% Input:
% varname  = name of variable (string)
% name	   = name of entry to load
% namelist = list of saved names
% var1mat  = matrix of saved variables (vectros in columns)
% var2mat  = 2nd matrix of saved variables (optional)
%
% Output:
% var1	   = variable to save (vector)
% var2	   = 2nd variable to save (optional)
%
% ver 3.3, Juan M. Rius, Jan 1997

function [var1,var2] = loadname(varname,name,namelist,var1mat,var2mat)

if ~length(namelist),
   warndlg([varname ': List empty'],'Warning:');
	var1 = []; if nargin > 4, var2 = []; end
	return;
end

inlist = findname(namelist,name);

if all(~inlist),	% If name does no exist in list
   text = [varname ': Name ''' name ''' not in list'];
   warndlg(text,'Warning:');
   
else,	% for necessary if name matches with more than one entry
	for i = find(inlist)',
		text = [varname ': Load ''' deblank(namelist(i,:)) ''' ?'];
      clk = questdlg(text,'Atention:','Yes','No','No');
      if strcmp(clk,'Yes'),
			var1 = var1mat(:,i);
			if nargin > 4, var2 = var2mat(:,i);
			end
		end
	end
end;