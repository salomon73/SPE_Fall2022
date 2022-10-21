	function r_beam = r_beam_Cryogenic_170(magnet,I)
%
% function r_beam = r_beam_Cryogenic_170(magnet,I)
%
%	Ex:		>> r_beam = r_beam_Cryogenic_170('cryogenic',[57.712  , 68.596 ,  107.780  , 107.780])
%
% Returns an estimate of the electron beam radius in the cavity of the GT140 gyrotron
%%
% Version:	    1.0		J.-P. Hogge		Feb. 27, 2008
%               1.1     J.-P. Hogge     Nov. 23, 2011   : Refurbished first prototype parameters
%               2.0     J.-P. Hogge     Aug. 10th, 2017 : Estimation based on elliptical integrals only (a bit more complicated, but more precise)
%               2.1     J.-P. Hogge     Jun. 7th, 2018  : corrected minor bug (so far only cryogenic magnet was implemented.
%
%   
%   magnet	  :  'cryogenic' or any other 170GHz magnet
%
%
% Procedure:
% ----------
%
% First the azimuthal component aphi of the vector potential is computed  at the cathode.
% Then the axial component of the magnetic field in the vicinity of the axis, and
% in the cavity region is estimated.
% By applying flux conservation, the beam radius is estimated.
%
%****************************************************************************************

%   rcat	= 0.05318;	    % Cathode center radius [m]     1st 1MW Prototype, cold dimension (2017)   
%    rcat	= 0.0533693;	% Cathode center radius [m]     1st 1MW Prototype, warm dimension (2017)
%	zcat	= 0.136;	    % Cathode center location [m]   Refurbished Prototype (2011)
%	zcav	= 0.4995;	    % Cavity location [m]
%    
    GT170_1MW_Geometry

% Compute vector potential aphi at cathode location
	
	[Aphi]  = B_Ellip_Cryogenic_170('aphi',magnet,I,r_cat,z_cat);

%	flux    = 2 * pi * rcat * Aphi;

% Compute axial magnetic field and its derivatives in cavity

%	[Bzcav,dBzcav,d2Bzcav] = B_on_axis_170('asg_modified_201206',I,zcav);

    Bzcav   = B_Ellip_Cryogenic_170('bz',     magnet,I,0,z_cav);
    dBzcav  = B_Ellip_Cryogenic_170('dbzdz',  magnet,I,0,z_cav);
    d2Bzcav = B_Ellip_Cryogenic_170('d2bzdz2',magnet,I,0,z_cav);
    
%
% Apply flux conservation to find beam radius in cavity 
% The exact flux at cathode (r_cat * Aphi) is equalled to the flux
% estimated with paraxial expansion.
%
	a 		    = 1/16 * d2Bzcav;
	b 		    = - Bzcav/2;
	c 		    = r_cat * Aphi;
	
	rbsq_plus	= 1./(2*a) * (-b + sqrt(b^2 - 4*a*c));  
	rbsq_minus	= 1./(2*a) * (-b - sqrt(b^2 - 4*a*c));
	
	r_beam		= sqrt(rbsq_minus);
    
%
% Estimation based on elliptical integrals only. 
% The following equation  is solved for r_beam:
% r_cat * Aphi(r_cat,z_cat) = r_beam * Aphi(r_beam,z_cav)
%
    raphi_cat   = r_cat * Aphi;
   
    myfun       = @(r,z_cav,I,raphi_cat) r*B_Ellip_Cryogenic_170('aphi',magnet,I,r,z_cav)-raphi_cat;  

    fun         = @(r) myfun(r,z_cav,I,raphi_cat);
    
    r_beam      = fzero(fun, 0.009);
    
    
    
    
    
