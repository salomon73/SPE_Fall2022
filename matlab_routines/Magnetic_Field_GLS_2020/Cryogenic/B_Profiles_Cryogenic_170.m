	function  B_Profiles_170(magnet,I)
%	
%    Call:      B_Profiles_170(magnet,I)
%    Example:   >>  I = [9.7000    7.7000   96.2499   91.0124   87.4723]
%               >>  B_Profiles_170('asg_modified_201206',I);
%
%               >>  I = [61.56    66.45   113.43   108.09]/3;
%               >>  B_Profiles_170('cryogenic',I);
%
%
% Returns magnetic field and beam parameters associated to the GT170 GHz magnet
%
%	Version :	1.0		J.-P. Hogge		Feb. 10, 2005
%
%   magnet	: 'fzk_165ghz_oi'
%	I       : Array of currents, as read on PS
%	Dimensions are in [m]
	

%	z_shim       	= 0.000;
%	z_cat_design 	= 0.08;
%	z_cav_design 	= 0.434;
%	z_cat 			= z_cat_design ;
%	z_cav 			= z_cav_design ;
%	r_cat			= 0.05;

%   z_shim       	= 0.000;
% 	z_cat_design 	= 9999;
% 	z_cav_design 	= 0.394;
% 	z_cat 			= z_cat_design ;
% 	z_cav 			= z_cav_design ;
% 	r_cat			= 0.05;

% Cryogenic magnet for the EU 170GHz 1MW prototype
    z_shim       	= 0.000;
 	z_cat_design 	= 0.136;
 	z_cav_design 	= 0.499665;
 	z_cat 			= z_cat_design ;
 	z_cav 			= z_cav_design ;
 	r_cat			= 0.05318;


%	z				= linspace(0.00,0.900,91)-0.00012;             % =.00012 is the difference between reading and actual
                                                                    %  position of the probe (due to the black ring curvature)
                                                                    % JPH: 2012_01_12 
                                                                    
	z				= linspace(0.00,0.900,901)-0.000;             % 

%	z_interest		= [z_cat, z_cav, 0.228,0.296,0.358,0.512,0.572];				% A few points with particular interest
    z_interest      = [];
    z				= cat(2,z,z_interest);									% Add those locations to the linspace
	z				= sort(z);												% And finally sort
	r				= 0 * z;
	
%	'bz','br','dbzbz','dbzdr','dbrdz','dbrdr', 'd2bzdz2' or 'aphi'

	Bz				=  B_Ellip_Cryogenic_170('bz',     magnet,I,r,z);
	dBzdz			=  B_Ellip_Cryogenic_170('dbzdz',  magnet,I,r,z);
	d2Bzdz			=  B_Ellip_Cryogenic_170('d2bzdz2',magnet,I,r,z);
%
	r				=  r+0.025;		
	Bz_25			=  B_Ellip_Cryogenic_170('bz',     magnet,I,r,z);
	Br_25			=  B_Ellip_Cryogenic_170('br',     magnet,I,r,z);
	dBzdz_25		=  B_Ellip_Cryogenic_170('dbzdz',  magnet,I,r,z);
	dBzdr_25		=  B_Ellip_Cryogenic_170('dbzdr',  magnet,I,r,z);

%	Br				=  B_Ellip_140('br',magnet,Igun,Icav,r,z)
%	dBzdz			=  B_Ellip_140('dbzdz',magnet,Igun,Icav,r,z)
%	dBzdr			=  B_Ellip_140('dbzdr',magnet,Igun,Icav,r,z)
%	dBrdz			=  B_Ellip_140('dbrdz',magnet,Igun,Icav,r,z)
%	dBrdr			=  B_Ellip_140('dbrdz',magnet,Igun,Icav,r,z)
%	d2Bzdz			=  B_Ellip_140('d2bzdz2',magnet,Igun,Icav,r,z)
%	aphi			=  B_Ellip_140('aphi',magnet,Igun,Icav,r,z)
	
%	r_beam	=  r_beam_140(magnet,Igun,Icav);

	disp(' ')	
%	disp(sprintf('Cathode magnetic field on axis (r=0) :  %1.5f [T]',Bz(1)))
%	disp(sprintf('Cathode magnetic field at r=r_cat    :  %1.5f [T]',Bz(2)))
%	disp(sprintf('Cavity magnetic field                :  %1.5f [T]',Bz(3)))
%	disp(sprintf('Magnetic field compression  (r=0)    : %2.5f',Bz(3)/Bz(1)))
%	disp(sprintf('Magnetic field compression  (r=rcat) : %2.5f',Bz(3)/Bz(2)))
%	disp(sprintf('Beam radius in cavity                : %2.5f [mm]' ,r_beam*1000))
%	disp(sprintf('Beam Compression (r_cat/r_beam)^2    : %1.5f' ,(r_cat/r_beam)^2))
	disp(' ')
	
	
		         
	disp(' ')
	disp(sprintf('GT170 GHz magnetic field and derivatives for I1,..,I5=%0.5g %0.5g %0.5g %0.5g  [A]',I(1),I(2),I(3),I(4)))
	disp('-------------------------------------------------------------------------------')
	disp(' ')
	disp(sprintf('Cathode is located in z_cat=%0.5g [m]',z_cat))
	disp(sprintf('Cavity  is located in z_cav=%0.5g [m]',z_cav))
	disp(' ')	
	disp('z [m]   Bz [T]   dBz/dz [T/m]    d2Bz/dz2 [T/m2]   Bz(r=25mm)   Br(r=25mm)   dBz/dz(r=25mm)  dBz/dr(r=25mm)  ')
	disp('-----   ------   ------------    ---------------   ----------   ----------   --------------  --------------  ')

	for i=1:length(z)
		disp(sprintf('%8.5f   %8.5f    %9.5f     %10.5f     %8.5f      %8.5f     %9.5f     %9.5f    ',...
		              z(i),   Bz(i),  dBzdz(i),  d2Bzdz(i), Bz_25(i), Br_25(i), dBzdz_25(i), dBzdr_25(i)   ))
	end

	disp(' ')	


%
%--------------------------------------------------------------------------------------------------
%
%	
	figure(1)

	subplot(3,1,1)
		plot(z,Bz);

		title(sprintf('GT170: Magnetic field with I1,..,I4=%0.5g %0.5g %0.5g %0.5g  [A]\n',I(1),I(2),I(3),I(4)))
		xlabel('Position z [m]')
		ylabel('B_z[T]')
		grid on	

	subplot(3,1,2)
		plot(z,dBzdz);

		title(sprintf('GT170: Magnetic field derivative with I1,..,I4=%0.5g %0.5g %0.5g %0.5g  [A]\n',I(1),I(2),I(3),I(4)))
		xlabel('Position z [m]')
		ylabel('dB_z/dz[T/m]')
		grid on	

	subplot(3,1,3)
		plot(z,d2Bzdz);

		title(sprintf('GT170 Magnetic field second derivative with I1,..,I4=%0.5g %0.5g %0.5g %0.5g  [A]\n',I(1),I(2),I(3),I(4)))
		xlabel('Position z [m]')
		ylabel('d2B_z/dz2[T/m2]')	
		grid on

%
%--------------------------------------------------------------------------------------------------
%

	filename = sprintf('b170profile_%s.txt',magnet)

	fid = fopen(filename,'w')
	
	fprintf(fid,'GT170 GHz magnetic field and derivatives for I1,..,I5=%0.5g %0.5g %0.5g %0.5g  [A]\n',I(1),I(2),I(3),I(4));
	fprintf(fid,'------------------------------------------------------------- \n');
	fprintf(' \n');
	fprintf(fid,'Cathode is located in z_cat=%0.5g [m] \n',z_cat);
	fprintf(fid,'Cavity  is located in z_cav=%0.5g [m] \n',z_cav);
	fprintf(' \n');
	fprintf(fid,' z [m]   Bz [T]   dBz/dz [T/m]    d2Bz/dz2 [T/m2]   Bz(r=25mm)   Br(r=25mm)   dBz/dz(r=25mm)  dBz/dr(r=25mm)   \n');
	fprintf(fid,' -----   ------   ------------    ---------------   ----------   ----------   --------------  --------------   \n');

	for i=1:length(z)
		fprintf(fid,'%8.5f   %8.5f    %9.5f     %10.5f     %8.5f      %8.5f     %9.5f     %9.5f   \n',...
		             z(i),   Bz(i),  dBzdz(i),  d2Bzdz(i), Bz_25(i), Br_25(i), dBzdz_25(i), dBzdr_25(i)    );
	end
	
	fclose(fid)
	
