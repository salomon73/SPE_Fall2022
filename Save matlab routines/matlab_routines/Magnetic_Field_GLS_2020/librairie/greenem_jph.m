function b = greenem_jph(mode,r1,z1,r2,z2)
%
%	Call:	b = greenem_jph(mode,r1,z1,r2,z2)
%
% Green's functions
%
% The Green's functions are estimated at the position of the system 1 for a unit current flowing in system 2, 
% for example BR_1_2 is the radial magnetic field at position 1 produced by 1A at position 2. 
% In the formulae given bellow, system 1 is at position (r1,z1) and system 2 at position (r2,z2).
%
%
%	Version:	1.0		J.-P. Hogge		Oct. 4, 2001  (based on an original routine written by J.-M. Moret)
%				1.1		J.-P. Hogge		Jul. 10 2002  Added aphi as an option for the mode
%				1.2		J.-P. Hogge		Nov. 5  2002  Corrected mistake in dbzdr calculation
%
%   mode    : 'bz','br','dbzbz','dbzdr','dbrdz','dbrdr', 'd2bzdz2' or 'aphi'
%   r1      : 1D array of radial       locations where field is computed
%   z1      : 1D array of longitudinal locations where field is computed
%   r2      : 1D array of radii        of the subcoils producing the field
%   z2      : 1D array of z-locations  of the subcoils producing the field
%
%
%
%       Definitions
%       -----------
%       h   = z2 - z1
%       d^2 = (r2 - r1)^2 + h^2
%       u^2 = (r2 + r1)^2 + h^2
%       k^2 = 4 * r1 * r2 / u^2
%       v^2 = r2^2 + r1^2 + h^2
%       w^2 = r2^2 - r1^2 - h^2
%
%       E(k^2) and K(k^2) and the complete elliptic integrals of the first and second kind respectively 
%
%		Magnetic field
%       --------------
%       br = -(1/2/pi/r1) * (dmut/dh)  = (mu0/2/pi) * (h/r1/u) * ((v^2 / d^2) * E(k^2) - K(k^2))
%       bz =  (1/2/pi/r1) * (dmut/dr1) = (mu0/2/pi) * (1/u)    * ((w^2 / d^2) * E(k^2) + K(k^2))
%       bt =  br * cos(ang) + bz * sin(ang)
%       bn =  bz * cos(ang) - br * sin(ang)
%
%       Magnetic field derivatives
%       --------------------------
%       d(bz)/dr1 = (mu0/2/pi) * (1/r1/d^2/u^3) * ((v^2*d^2*u^2 - h^2*(d^4+k^2*u^4))/d^2 * E(k^2) +              (h^2*v^2-u^2*d^2) * K(k^2))
%       d(bz)/dz1 = (mu0/2/pi) * (   h/d^2/u^3) * (           (-4*v^2*w^2-3*u^2*d^2)/d^2 * E(k^2) +                            w^2 * K(k^2))
%       d(br)/dr1 = -d(bz)/dz1 - br/r1
%       Here is a small change w respect to J-M Moret's original routine
%       d(br)/dz1 = (mu0/2/pi) * ( h^2/d^2/u/r1) * ( (-u^2/d^2 - d^2/u^2 + v^2/h^2 +1)   * E(k^2) +           (-d^2/h^2 - v^2/u^2) * K(k^2))




	mu0	= 4 * pi * 1.0e-7;


% Since we want to know the fields at locations (r1,z1) produced by sources
% at locations (r2,z2), we construct rectangular matrices of dimension
% (length(r1), length(r2)).

	if size(r1,1) > 1 | size(r2,1) > 1  , error('Input must be row vectors'); end
	
	R1	= repmat(r1', 1         , length(r2) );		% R1 is a matrix length(r1)xlength(r2)
	Z1	= repmat(z1', 1         , length(r2) );		% Z1 is a matrix length(r1)xlength(r2)
	R2	= repmat(r2 , length(r1), 1          );		% R2 is a matrix length(r1)xlength(r2)
	Z2	= repmat(z2 , length(r1), 1          );		% Z2 is a matrix length(r1)xlength(r2)
		
	
	h	= Z1 -Z2;
	d	= sqrt( (R2 - R1).^2 + h.^2 );
	u	= sqrt( (R2 + R1).^2 + h.^2 );
	k	= sqrt(  4*R1.*R2 ./ u.^2 );
	v	= sqrt(  R2.^2 + R1.^2 + h.^2);
	w	= sqrt(  R2.^2 - R1.^2 - h.^2);
	
	
	[K,E] = ellipke(k.^2);
	
	
	if(~iscell(mode))
        mode={mode};
    end
    b=zeros([size(k),length(mode)]);
        for i=1:length(mode)
        
        switch lower(mode{i})

 
 	     case 'bz', 

			b(:,:,i)    =  (mu0/2/pi) * (1./u)      .* (  (w.^2./d.^2).*E + K  );
			 
	     case 'br', 
	 
			b(:,:,i)	=   (mu0/2/pi) * (h./R1./u)  .* (  (v.^2./d.^2).*E - K  );
	 			 
	  	 case 'dbzdz',
	 
			b(:,:,i)	=   (mu0/2/pi) * (h./d.^2./u.^3)  .* (  (-4*v.^2.*w.^2-3*u.^2.*d.^2)./d.^2.*E   +   w.^2 .* K  );	 
			 
		 case 'dbzdr', 
	 
%JPH 5_11_2002			b	=   (mu0/2/pi) * (1./R2./d.^2./u.^3)  .* ( (v.^2.*d.^2.*u.^2 - h.^2.*(d.^4+k.^2.*u.^4))./d.^2.*E   +   (h.^2.*v.^2 - u.^2.*d.^2) .* K  );	 
			b(:,:,i)	=   (mu0/2/pi) * (1./R1./d.^2./u.^3)  .* ( (v.^2.*d.^2.*u.^2 - h.^2.*(d.^4+k.^2.*u.^4))./d.^2.*E   +   (h.^2.*v.^2 - u.^2.*d.^2) .* K  );	 
			 
		 case 'dbrdr', 
	 
%        d(br)/dr1 = -d(bz)/dz1 - br/r1

			b(:,:,i)   = - (mu0/2/pi) * (h./d.^2./u.^3)  .* (  (-4*v.^2.*w.^2-3*u.^2.*d.^2)./d.^2.*E   +   w.^2 .* K  ) ...
			      - (mu0/2/pi) * (h./R1./u)  .* (  (v.^2./d.^2).*E - K  ) ./R1;
			 
		 case 'dbrdz', 
	 
			b(:,:,i)	=   (mu0/2/pi) * (h.^2./d.^2./R1./u)  .* ( ( -u.^2./d.^2 - d.^2./u.^2 + v.^2./h.^2 +1).*E          +  (-d.^2./h.^2 + v.^2./u.^2) .* K  );	 
			 
	     case 'd2bzdz2', 
	 
			b(:,:,i) 	=	(mu0/2/pi) * (1./d.^4./u.^3) .* (    (           ( -d.^2.*u.^2 + 4*h.^2.*( 2*R1.*R2 + d.^2)) + 8 * h.^2.*u.^2.*w.^2./d.^2   ...
														   + 2 *     ( -u.^2.*w.^2 + 4*h.^2.*( u.^2  + w.^2))    - 2 * d.^2.*( u.^2  + w.^2 )   ...
														   - h.^2 .* ( -8*d.^2     - (7 - 8*k.^2).*w.^2) )  .* E                                ...
													   - (              d.^2.*h.^2 + 4 * h.^2.*w.^2  + 											...
													   d.^2./u.^2 .* ( -u.^2.*w.^2 + 4*h.^2.*( u.^2  + w.^2)) ) .* K );
			 

		case 'aphi'											   
													   
			b(:,:,i)   =   (mu0/pi)   * sqrt(R2./(k.^2.*R1))  .* ( (1-k.^2/2).*K - E);									   
		otherwise, 
			disp('Unknown mode')
	  
		  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			  
		end
        end


return
