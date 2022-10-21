	function [M] = B_Matrix_Cryogenic_170(magnet)
%	
%    Call:        [M] = B_Matrix_Cryogenic_170(Magnet)
%    Example:  >> [M] = B_Matrix_Cryogenic_170('asg_modified_201206');
%
%	Version :	1.0		J.-P. Hogge		Mar. 27, 2008
%               1.1     J.-P. Hogge     Nov. 23, 2011, Changet Cathode parameters to refurbished first prototype
%
%
% Magnetic field matrix associated with GT170 GHz magnet. 
% This scripts determines quantities M(1,1)...M(4,4) such that
%   
%
%  B_gun            M(1,1)  M(1,2)  M(1,3)  M(1,4)   M(1,5)   I1
%              =                          *         
%  B_cavity         M(2,1)  M(2,2)  M(2,3)  M(2,4)   M(2,5)   I2
%
% dBz/dz_gun        M(3,1)  M(3,2)  M(3,3)  M(3,4)   M(3,5)   I3
%
% dBz/dz_cav        M(4,1)  M(4,2)  M(4,3)  M(4,4)   M(4,5)   I4
%
%
% with currents I1, -I2, -I3, I4  flowing in the 4 power supplies
%
	
% GT170 1MW first prototype (TH1509) geometry
% %
% 	z_shim       	= 0.00;
% 
%     z_cat_design 	= 0.136;                    % emitter radius AT CENTER
%     z_cav_design 	= 0.4995;                   % cavity center position
%     z_launcher_d    = 0.78848;                  % launcher edge 
%     
% 	z_cat 			= z_cat_design + z_shim;
% 	z_cav 			= z_cav_design + z_shim;
%     z_launcher      = z_launcher_d + z_shim;
% 
%     r_cat			= 0.053369;                 % emitter radius AT CENTER, WARM CONDITION
%     r_launcher      = 0.02068;                  % launcher radius
    
%  Include GT170 Geometry
    
    GT170_1MW_Geometry
%
%
 	r		= [0     0       ] % 
   	z		= [z_cat z_cav   ]


% this script is specially dedicated to the 170GHz Cryogenic magnet...
%%%% JHE 2018_11_27   magnet = 'cryogenic' 
    disp(sprintf('Aimant  : %s',magnet))
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[1 0 0 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[1 0 0 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[1 0 0 0],r,z);    

    column   = 1;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=DB(2);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 1 0 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 1 0 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 1 0 0],r,z);
    
    column   = 2;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=DB(2);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 0 1 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 0 1 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 0 1 0],r,z);

    column   = 3;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=DB(2);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 0 0 1],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 0 0 1],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 0 0 1],r,z);

    column   = 4;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=DB(2);
    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
%     [B] 	= B_Ellip_170('bz'     ,magnet,[0 0 0 0 1],r,z);
% 	[DB] 	= B_Ellip_170('dbzdz'  ,magnet,[0 0 0 0 1],r,z);
% 	[D2B] 	= B_Ellip_170('d2bzdz2',magnet,[0 0 0 0 1],r,z);
%     
%     column   = 5;
%     M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=DB(2); M(5,column)=D2B(1);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp(sprintf('  Bz_gun            %f     %f     %f     %f            I1 ', M(1,1), M(1,2), M(1,3), M(1,4)))
    disp(sprintf('  Bz_cav            %f     %f     %f     %f            I2 ', M(2,1), M(2,2), M(2,3), M(2,4)))
    disp(sprintf('  dBzdz_gun         %f     %f     %f     %f            I3 ', M(3,1), M(3,2), M(3,3), M(3,4)))
    disp(sprintf('  dBzdz_cav         %f     %f     %f     %f            I4 ', M(4,1), M(4,2), M(4,3), M(4,4)))
    
	return
