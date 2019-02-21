% [namelist,var1mat,var2mat] = delname(varname,name,namelist,var1mat,var2mat)
% Delete results from list. If name matches with an already saved name,
% a dialog window prompts for confirmation
%
% Input:
% varname  = name of variable (string)
% name	   = name of entry to delete
% namelist = list of saved names
% var1mat  = matrix of saved variables (vectros in columns)
% var2mat  = 2nd matrix of saved variables (optional)
%
% Output:
% Updated variables
%
% ver 3.3, Juan M. Rius, Jan 1997

function [namelist,var1mat,var2mat] = delname(varname,name,namelist,var1mat,var2mat)

N = size(namelist,1);	% Number of entries in list

if N==0,	% List empty
   warndlg([varname ': List empty'],'Warning:');
 	return;
end

inlist = findname(namelist,name);
namenew = []; var1new = []; var2new = [];

% for loop necessary if more than one entry matches with name
% new variables are necessary because deleting a row or column changes
% the indices of the next rows and columns
for i = 1:N,
	word = deblank(namelist(i,:));
	if ~inlist(i),	% Do not delete this entry
		if ~length(namenew), namenew = word;
		else namenew = str2mat(namenew,word);
		end
		var1new = [var1new var1mat(:,i) ];
		if nargin > 4, var2new = [var2new var2mat(:,i) ];
		end
	else
		text = [varname ': Delete ''' word ''' ?'];
		clk = questdlg(text,'Atention:','Yes','No','Cancel','No');

      switch clk,
      case 'No',
			if ~length(namenew), namenew = word;
			else namenew = str2mat(namenew,word);
			end
			var1new = [var1new var1mat(:,i) ];
			if nargin > 4, var2new = [var2new var2mat(:,i) ];
         end
      case 'Cancel',
         return
		end
	end
end

namelist = namenew; var1mat = var1new; if nargin>4, var2mat = var2new; end



