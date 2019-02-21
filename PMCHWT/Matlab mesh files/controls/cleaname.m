% [namelist,var1mat,var2mat] = cleaname(varname,namelist,var1mat,var2mat)
% Delete all results from list. A dialog window prompts for confirmation
%
% Input:
% varname  = name of variable (string)
% namelist = list of saved names
% var1mat  = matrix of saved variables (vectros in columns)
% var2mat  = 2nd matrix of saved variables (optional)
%
% Output:
% Updated variables (all [])
%
% ver 3.3, Juan M. Rius, Jan 1997

function [namelist,var1mat,var2mat] = cleaname(varname,namelist,var1mat,var2mat)

if ~length(namelist),	% List empty
   warndlg([varname ': List already empty'],'Warning:');
	namelist = []; var1mat = [];
	if nargin > 3, var2mat = [];
	end 
	return;
end

clk = questdlg([varname ': Clear all saved results?'],'Atention:','Yes','No','No');

if strcmp(clk,'Yes'),
	namelist = []; var1mat = [];
	if nargin > 3, var2mat = [];
	end 
end

