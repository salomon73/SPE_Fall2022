	function [M,gamma,delta] = B_Mat88_Cryogenic_170(magnet)
%	
%    Call:        [M,gamma,delta] = B_Mat88_Cryogenic_170(Magnet)
%    Example:  >> [M,gamma,delta] = B_Mat88_Cryogenic_170('toto');
%
%	Version :	1.0		J.-P. Hogge		Aug. 3rd, 2017
%                       largely based on B_Matrix_170.m
%
% Magnetic field matrix associated with GT170 GHz Cryogenic magnet. 
% This a set of coefficient such that their scalar product with the
% currents provided by the power supplies should be < 0 to avoid
% interception. The gamma parameters correspond to the case where magnetici
% flux conservation is used, whereas the delty parameters are based on the
% vector potential A_phi.
%
%   
%
%  B_gun            M(1,1)  M(1,2)  M(1,3)  M(1,4)   M(1,5)   I1
%              =                          *         
%  B_cavity         M(2,1)  M(2,2)  M(2,3)  M(2,4)   M(2,5)   I2
%
% dBz/dz_gun        M(3,1)  M(3,2)  M(3,3)  M(3,4)   M(3,5)   I3
%
%  B_launcher       M(4,1)  M(4,2)  M(4,3)  M(4,4)   M(4,5)   I4
%
%
% with currents I1, I2, -I3, I4  flowing in the 4 power supplies
%
	
% GT170 1MW first prototype (TH1509) geometry
%
	z_shim       	= 0.00;

    z_cat_design 	= 0.1335;                   % corresponds to the LARGEST emitter radius,
                                                % not the average one.
    z_cav_design 	= 0.4995;                   % cavity center position
    z_launcher_d    = 0.78848;                  % launcher edge 
    
	z_cat 			= z_cat_design + z_shim;
	z_cav 			= z_cav_design + z_shim;
    z_launcher      = z_launcher_d + z_shim;

    r_cat			= 0.05419;                  % largest emitter radius
    r_launcher      = 0.02068;                  % launcher radius
    
    
	

%	r		= [0     0       0      0]          % If one decides to use on-axis values
 	r		= [0     0       0      r_launcher] % If one decides to use more stringent values
   
    rphi    = [r_cat 0       r_cat  r_launcher] % check with vector potential
	z		= [z_cat z_cav   z_cat  z_launcher]


% this script is specially dedicated to the 170GHz Cryogenic magnet...
    magnet = 'cryogenic' ;
    disp(sprintf('Aimant  : %s',magnet))
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[1 0 0 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[1 0 0 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[1 0 0 0],r,z);    
    [Aphi]  = B_Ellip_Cryogenic_170('aphi',   magnet,[1 0 0 0],rphi,z);

    column   = 1;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=B(4);
    
    aphig(column) = Aphi(1);
    aphil(column) = Aphi(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 1 0 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 1 0 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 1 0 0],r,z);
    [Aphi]  = B_Ellip_Cryogenic_170('aphi',   magnet,[0 1 0 0],rphi,z);

    
    column   = 2;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=B(4);
    
    aphig(column) = Aphi(1);
    aphil(column) = Aphi(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 0 1 0],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 0 1 0],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 0 1 0],r,z);
    [Aphi]  = B_Ellip_Cryogenic_170('aphi',   magnet,[0 0 1 0],rphi,z);

    column   = 3;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=B(4);
    
    aphig(column) = Aphi(1);
    aphil(column) = Aphi(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    [B] 	= B_Ellip_Cryogenic_170('bz'     ,magnet,[0 0 0 1],r,z);
	[DB] 	= B_Ellip_Cryogenic_170('dbzdz'  ,magnet,[0 0 0 1],r,z);
	[D2B] 	= B_Ellip_Cryogenic_170('d2bzdz2',magnet,[0 0 0 1],r,z);
    [Aphi]  = B_Ellip_Cryogenic_170('aphi',   magnet,[0 0 0 1],rphi,z);

    column   = 4;
    M(1,column)  = B(1) ; M(2,column) = B(2); M(3,column)=DB(1); M(4,column)=B(4);
    
    aphig(column) = Aphi(1);
    aphil(column) = Aphi(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp(sprintf('  Bz_gun         %f     %f     %f     %f            I1 ', M(1,1), M(1,2), M(1,3), M(1,4)))
    disp(sprintf('  Bz_cav         %f     %f     %f     %f            I2 ', M(2,1), M(2,2), M(2,3), M(2,4)))
    disp(sprintf('  dBzdz_gun      %f     %f     %f     %f            I3 ', M(3,1), M(3,2), M(3,3), M(3,4)))
    disp(sprintf('  Bz_launcher    %f     %f     %f     %f            I4 ', M(4,1), M(4,2), M(4,3), M(4,4)))
    
    gamma = [0 0 0 0];
    delta = [0 0 0 0];

    for i=1:4
        gamma(i) = M(1,i) - (r_launcher/r_cat)^2 * M(4,i);
        delta(i) = aphig(i) - (r_launcher/r_cat) * aphil(i)
    end
    disp(' ')
    format long
    disp(sprintf('  gamma         %f     %f     %f     %f     ',gamma(1), gamma(2), gamma(3), gamma(4)))
    disp(sprintf('  delta         %f     %f     %f     %f     ',delta(1), delta(2), delta(3), delta(4)))

 
	return
