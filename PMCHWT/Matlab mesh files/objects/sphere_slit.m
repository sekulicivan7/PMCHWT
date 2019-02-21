 % composite object of two cubes or a cube and a pyramid one on the another sharing the
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lam	= 1;
    k = 2*pi/lam;
   eta= 120*pi; 	
    mu0=4*pi*1e-7; 	
   eps0=8.8541878176*1e-12; 
   epsr1=1;
   mur1=1;
    n_S_f = 9;	
	n_S_s = 9;
 	n_V_f = 11;	
 	n_V_s = 11;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
       th_i	= 180;	
	   ph_i	= 0;	
	   rot_EH  = 270; 
       
    ang_e1 = 0;
	ang_e2 = 180;	
    M_e = 60;	
	ang_h1 = 0;	
	ang_h2 = 180;	
	M_h = 60;

 geom	= 'sphere_Ivan';	
       param1 = [0.1 192 1];
        param=param1;
       fact_vol=1/100; 
       run('run_geom_Ivan')
       obj1=obj;


       r1 = obj1.vertex( : , obj1.topol(1,:) );
       r2 = obj1.vertex( : , obj1.topol(2,:) );
       r3 = obj1.vertex( : , obj1.topol(3,:) );
       
       k1=find(r1(3,:)>=0);
       k2=find(r2(3,:)>=0);
       k3=find(r3(3,:)>=0);

  Z11=intersect(k1,k2);
       UP=intersect(Z11,k3);  % ovi trokuti dijele dvije regije
       
       topolUP=obj1.topol(:,UP);
       vertexUP=obj1.vertex;
       triangleUP=obj1.trian(:,UP);
       unUP=obj1.un(:,UP);
       dsUP=obj1.ds(UP);
       centUP=obj1.cent(:,UP);

 objUP = struct('vertex',[vertexUP],'topol',[topolUP],'trian',[triangleUP],'un',[unUP],'ds',[dsUP],'cent',[centUP]);


       obj_nc_inpVUP = get_obj_nc_inprismV_Ivan( objUP , fact_vol );
       obj_nc_inpV2UP = get_obj_nc_inprismV_Ivan2( objUP , fact_vol );



%trazenje po drugoj regiji
 geom	= 'sphere_Ivan';	
       param2 = [0.1 192 1];
        param=param2;
       run('run_geom_Ivan')
       obj2=obj;


       r1 = obj2.vertex( : , obj2.topol(1,:) );
       r2 = obj2.vertex( : , obj2.topol(2,:) );
       r3 = obj2.vertex( : , obj2.topol(3,:) );
       
       k1=find(r1(3,:)<=0);
       k2=find(r2(3,:)<=0);
       k3=find(r3(3,:)<=0);
       
       Z12=intersect(k1,k2);
       BOT=intersect(Z12,k3);%ovi trokuti dijele dvije regije
       
       topolBOT=obj2.topol(:,BOT);
       vertexBOT=obj2.vertex;
       triangleBOT=obj2.trian(:,BOT);
       unBOT=obj2.un(:,BOT);
       dsBOT=obj2.ds(BOT);
       centBOT=obj2.cent(:,BOT);

       objBOT = struct('vertex',[vertexBOT],'topol',[topolBOT],'trian',[triangleBOT],'un',[unBOT],'ds',[dsBOT],'cent',[centBOT]);


       obj_nc_inpVBOT = get_obj_nc_inprismV_Ivan( objBOT , fact_vol );
       obj_nc_inpV2BOT = get_obj_nc_inprismV_Ivan2( objBOT , fact_vol );
       