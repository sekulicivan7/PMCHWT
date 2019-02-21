% Functions for creating uicontrols that modify variables
%
% Version 3.2, October 1996, Juan M.Rius
% Version 3.3, January 1997, Juan M. Rius
%
% a) This functions allow creation of uicontrols that modify variables
%    in global space or current window UserData
% b) Callbacks do not set any variable in global space!!
% c) v3.3, init_win can be called from callbacks, to allow creation of dlgwin  
%
% Main functions for creating figure with uicontrols:
% init_win	- Initialize figure for uicontrols
%		  Get figure information necessary for uicontrol functions
% grid_win	- Draw lines-characters grid in current window
% inputbox	- Create input box in current window
% pushbutt	- Create push button in current window
% radiobut	- Create radio button array in current window
% checkbox	- Create check box in current window
% popupbox	- Create popup box in current window
% stattext	- Create static text in current window
% slider	- Create slider in current window
% frame		- Create frame in current window
% ed_box	- Obsolete function. Equivalent to the above
% rectang	- Draws rectangle in new axis, in "character" units
% addcallb	- Append strings to callback of uicontrol
% dlgwin	- Creates dialog window
%
% Useful functions:
% gud		- Get UserData
% sud		- Set UserData
% view_var	- View nicely the contents of variable with prefix text
%
% Functions used by the main functions:
% setradio	- Set radio button to 'on' and rest of radio array to 'off'
% updatvis	- Updates visibility of handles according to current value
%		  of all radiobuttons and checkboxes in current figure
% visibil	- Calls updatvis to update visibility of handles that depend on uih
%
% Database for storing results:
% savename	- Saves result in database that matches with name
% loadname	- Loads result from database that matches with name
% delname	- Deletes result from database that matches with name
% cleaname	- Clears all results in database
% findname	- Finds result in database that matches with name
