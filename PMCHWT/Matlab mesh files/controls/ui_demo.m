% RUN.M
% Interface for programs running with 2-D PEC surfaces.
%
% IE-MEI, version 3.0, Juan M. Rius, sept. 1996

figure(1); set(gcf, 'Visible', 'Off'); 	clf;
init_win(86,21.5,'IE-MEI for 2-D PEC cylinders');

% Initialize non-existing variables

if ~exist('geom'),	geom = '';	end;
if ~exist('param'),	param=[];	end;
if ~exist('N'),		N=[];		end;
if ~exist('lam'),	lam=[];		end;
if ~exist('TM'),	TM=1;		end;
if ~exist('ang_i'),	ang_i=[];	end;
if ~exist('ang_1'),	ang_1=[];	end;
if ~exist('ang_2'),	ang_2=[];	end;
if ~exist('M'),		M=[];		end;
if ~exist('out_par'),	out_par=zeros(1,7); end;
if ~exist('save_name'), save_name='';	end;
if ~exist('strsave'),	strsave=[];	end
if ~exist('Jsave'),	Jsave=[];	end
if ~exist('RCSmsave'),	RCSmsave=[];	end
if ~exist('RCSbsave'),	RCSbsave=[];	end
if ~exist('LA'),	LA=[];		end;
if ~exist('LB'),	LB=[];		end;
if ~exist('LEn'),	LEn=[];		end;
if ~exist('tEn'),	tEn=[];		end;
if ~exist('thres'),	thres=[];	end;
if ~exist('P'),		P=[];		end;
if ~exist('use_A'),	use_A=1;	end; 
if ~exist('use_B'),	use_B=1;	end;
if ~exist('ab_t'),	ab_t=0;		end;
if ~exist('met_t'),	met_t=0;	end;
if ~exist('comp_t'),	comp_t=0;	end;
if ~exist('intype'),	intype=1;	end;
if ~exist('F'),		F=1;		end;
if ~exist('fpedge'),	fpedge=4;	end;

% Object group
stattext(2,1,'Object geometry');	% Group title
frame(2, 6.85, 40, 5.7);

objlist=str2mat('airfoil','cavity','circular','eliptic','ogival','oval','naca0012');
objlist=str2mat(objlist,'rectang','scavity','semi','strip','strip_h');
popupbox(5,3.55,'geom',objlist,'... or select object here: ',24,10);

inputbox(5,2.25,'geom',geom,  'Object boundary = ',24,10);
inputbox(5,4.5, 'param',param,'Geometry parameters =',19,15);
inputbox(5,5.5, 'N',N,        'Number of basis functions = ',28,6);
inputbox(5,6.5, 'lam',lam,    'Wavelength = ',28,6);

% Results group
stattext(2,8,'Current and RCS');				% Group title
frame(2, 12.85, 40, 4.7);
inputbox(5,9.5, 'ang_i',ang_i,'Incident angle = ',28,6);
inputbox(5,10.5,'ang_1',ang_1,'First observation angle = ',28,6);
inputbox(5,11.5,'ang_2',ang_2,'Last observation angle = ',28,6);
inputbox(5,12.5,'M',M,        'Number of observation angles = ',28,6);
radiobut([5;12],[9.5; 9.5],'TM',[TM;1;0],str2mat('TM','TE'));

% Output group
stattext(2,14,'Output');			% Group title
frame(2, 20.85, 40, 6.7);
checkbox(4, 15.5,'out_par(1)',out_par(1),'Boundary ');
checkbox(4, 16.5,'out_par(2)',out_par(2),'abs (J)');
checkbox(18,16.5,'out_par(3)',out_par(3),'angle (J)');
checkbox(4, 17.5,'out_par(4)',out_par(4),'Mono. RCS ');
checkbox(18,17.5,'out_par(5)',out_par(5),'Bista. RCS ');
checkbox(18,15.5,'out_par(6)',out_par(6),'Linear system');
checkbox(31,17.5,'out_par(7)',out_par(7),'Matrix data');
inputbox(4,18.75,'save_name',save_name,'Result legend =',13);

runsave = ['Jsave=[Jsave J]; RCSmsave=[RCSmsave RCSm]; RCSbsave=[RCSbsave RCSb];',...
	   'if strsave, strsave = str2mat(strsave,save_name); else strsave = save_name; end;'...
	  ];

rundel = ['[tmp1,tmp2]=size(strsave); tmp2=[]; ',...
	  'for tmp3=1:tmp1, if findstr(strsave(tmp3,:),save_name), tmp2=[tmp2 tmp3]; end; end; ',...
	  'tmp3=ones(1,tmp1); tmp3(tmp2)=zeros(size(tmp2)); tmp2=find(tmp3); ',...
	  'strsave=strsave(tmp2,:); ',...
	  'Jsave=Jsave(:,tmp2); RCSmsave=RCSmsave(:,tmp2); RCSbsave=RCSbsave(:,tmp2); ',...
	  'clear tmp1 tmp2 tmp3',...
	 ];

runclear = 'Jsave=[]; RCSmsave=[]; RCSbsave=[]; strsave=[];';

pushbutt(33,16   ,[],[],'output','Plot');
pushbutt( 5,20.25,[],[],runsave, 'Add');
pushbutt(15,20.25,[],[],rundel,  'Remove');
pushbutt(27,20.25,[],[],runclear,'Clear all');
	
% MEI group.
stattext(44,1,'IE-MEI coefficients');			% Group title
frame(44, 10, 40, 9);
stattext(49,2.25,'[A][Es] + [B][Hs] = [En][Js]',[0 0 0.76]);	% Equation
edLA =	inputbox(55,3.5,'LA',LA,'Bandwidth = ',10,3);
edLB =	inputbox(55,4.5,'LB',LB,'Bandwidth = ',10,3);
cheuseA = checkbox(46,3.5,'use_A',use_A,'Use A',[edLA 1]);
cheuseB = checkbox(46,4.5,'use_B',use_B,'Use B',[edLB 1]);

stamet = stattext(46,7,'Metrons:');
edmet = inputbox(70,7,'P',P,'number = ',8,4);
radmet = radiobut([60;53],[7.1; 7.1],'met_t',[met_t;0;1],str2mat('Harmonic','Delta'),[edmet 0]);

edthr = inputbox(65,9.5,'thres',thres,'Trunc threshold = ',13,4);
edtEn = inputbox(68,9.5,'tEn',tEn,    'Trunc band = ',10,4);
stause = stattext(46,8.5,'Use En:');
raduse = radiobut([54;61;68;76],[8.5; 8.5; 8.5; 8.5],'comp_t',[comp_t;0;1;2;3],str2mat('Zero','Full','Sparse','Banded'),...
		[edthr 2 edtEn 3]);

edLEn = inputbox(75,6,'LEn',LEn,'= ',2,4);
stattext(46,6,'Desired En:');
raddes = radiobut([56;68;62],[6.1; 6.1; 6.1],'ab_t',[ab_t;0;1;2],str2mat('Zero','Banded','Full'),...
		[edLEn 1 edmet stamet radmet edLA edLB edthr edtEn 0 1]);

% Interpolation group

stattext(44,11.75,'Frequency extrapolation of [A],[B]');			% Group title
frame(44, 16.5, 40, 4.5);

inputbox(46,13.5,'F',F,'Extrapolation factor = ',18,5);
inputbox(46,14.5,'fpedge',fpedge,'Fixed points at edges = ',18,5);

stattext(46,16,'Interpolation type');
radiobut([60;67;75],[16; 16; 16],'intype',[intype;0;1;2],str2mat('FFT','Linear','Spline'));

runmei30 = ['ko = 2*pi/(lam*F); No = N/F; zl = feval(geom, No, param);',...
	    '[A, B, En, err, cnd] = mei_v30(zl, LA, LB, LEn, tEn, thres, ko, out_par(6), TM, P, ab_t, comp_t, met_t, use_A, use_B);'...
	   ];

runinter = ['t1 = clock; disp(''Interpolation'');',...
	    'A = int_edge(A, LA, zl, [F intype fpedge]);',...
	    'B = int_edge(B, LB, zl, [F intype fpedge]);',...
	    't2 = clock; see_var(''Time'', etime(t2,t1));',...
	    'zl = feval(geom, N, param);'...
	  ];

runjrcs = ['[J, RCSm, RCSb] = j_rcs(zl, A, B, En, ang_i, ang_1, ang_2,',...
	   'M, 2*pi/lam, TM, out_par(4), out_par(5), comp_t);',...
	  ];

pushbutt(73, 4, [], [],runmei30,'Compute');
pushbutt(72,14, 10, [],runinter,'Extrapolate');
pushbutt(4, 11, [], [],runjrcs, 'Compute');

% Credits group
stattext(46,18.5,'Author:');			% Group title
frame(44, 20.85, 40, 3.5);
stattext(50,19.5,'Juan M. Rius',[0.76 0 0]); 		% Credits
stattext(50,20.5,'Version 3.0, sept 1996',[0 0 0.76]); % Credits
pushbutt(76, 19.5, [], [],'close(gcf)','Quit');

% Final arrangements
set(gcf,'Visible','On');

