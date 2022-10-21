	function  B_Params_Cryogenic_170(magnet,I)
%	
%    Call:         B_Params_Cryogenic_170(magnet,I)
%    Example:  >>  B_Params_Cryogenic_170('cryogenic',[61.56  66.45  113.43  108.09 ]);
% 
%
%
% Returns magnetic field and beam parameters associated to the GT140 GHz magnet
%
%	Version :	1.0		J.-P. Hogge		Oct. 08, 2002
%               1.1     J.-P. Hogge     Feb   3, 2005   Changed Call (Current in array I)
%
%   magnet	: 'oxf_1' only so far
%	I      	: Array of coil currents as read on power supply [A]
	

% 	z_shim       	= 0.00;
% %	z_cat_design 	= 0.1185; % First prototype (2008)
% 	z_cat_design 	= 0.1225; % Refurbished First prototype (2011)
% 
% 	z_cav_design 	= 0.4995;
% 	z_cat 			= z_cat_design + z_shim;
% 	z_cav 			= z_cav_design + z_shim;
% %	r_cat			= 0.057;    % First prototype (2008)
% 	r_cat			= 0.059;    % Refurbished First prototype (2011)
%                                 % Note Value given by C.Li?vin on 2011.11.23
%                                 % by e-mail:  Mesure ? froid #20?C : 58,705mm
%                                 %             diam?tre ? chaud : 58,705*(1+5.10-6*930)=58,978mm

%  Include GT170 Geometry
    
    GT170_1MW_Geometry
                               
                                
	r		= [    0 r_cat     0];
	z		= [z_cat z_cat z_cav];
	
%	'bz','br','dbzbz','dbzdr','dbrdz','dbrdr', 'd2bzdz2' or 'aphi'

	Bz		=  B_Ellip_Cryogenic_170('bz',      magnet,I,r,z);
	dBzdz	=  B_Ellip_Cryogenic_170('dbzdz',   magnet,I,r,z);
	d2Bzdz	=  B_Ellip_Cryogenic_170('d2bzdz2', magnet,I,r,z);
    Br		=  B_Ellip_Cryogenic_170('br',      magnet,I,r,z);
	aphi	=  B_Ellip_Cryogenic_170('aphi',    magnet,I,r,z);
	
	r_beam	=  r_beam_Cryogenic_170(magnet,I);

    disp(' ')
    disp(' ')	
	disp(sprintf('Cathode magnetic field on axis (r=0)            :  %2.5f [T]',Bz(1)))
	disp(sprintf('Cathode magnetic field at r=r_cat               :  %2.5f [T]',Bz(2)))
	disp(sprintf('Cavity magnetic field                           :  %2.5f [T]',Bz(3)))
    disp(' ')
	disp(sprintf('Cathode magnetic field derivative on axis (r=0) :  %2.5f [T/m]',dBzdz(1)))
	disp(sprintf('Cathode magnetic field derivative at r=r_cat    :  %2.5f [T/m]',dBzdz(2)))
	disp(sprintf('Cavity magnetic field derivative                :  %2.5f [T/m]',dBzdz(3)))
    disp(' ')
    disp(sprintf('Magnetic field compression  (r=0)               :  %2.5f',Bz(3)/Bz(1)))
	disp(sprintf('Magnetic field compression  (r=rcat)            :  %2.5f',Bz(3)/Bz(2)))
	disp(sprintf('Beam radius in cavity                           :  %2.5f [mm]' ,r_beam*1000))
	disp(sprintf('Beam Compression (r_cat/r_beam)^2               :  %2.5f' ,(r_cat/r_beam)^2))
	disp(' ')
	disp(sprintf('Magnetic angle at emitter [deg]                 :  %2.5f [deg]' ,atan2d(Br(2),Bz(2)) ))
    disp(' ')
	disp(' ')
	disp(' ')
	disp(' ')
	disp(' ')

	

	
		         
