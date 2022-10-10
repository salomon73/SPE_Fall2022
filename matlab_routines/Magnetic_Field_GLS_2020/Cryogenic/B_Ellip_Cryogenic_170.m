	function [B] = B_Ellip_Cryogenic_170(mode,magnet,I,r,z)
%	
%    Call:     [B] = B_Ellip_Cryogenic_170(mode,magnet,r,z)
%    Examples: 
%    1. Plot axis magnetic field profile:
%   
%              >> z = linspace(0,1,101);
%              >> r = 0*z + 0.0;
%              >> I = [61.56  66.45  113.43  108.09 ];      
%              >> [B] = B_Ellip_Cryogenic_170('bz','cryogenic',I,r,z);
%              >> figure(1); plot(z,B)
%
%    2. Plot field lines (r* aphi = cst)
%
%              >> z = linspace(0,1,51);
%              >> r = linspace(0,0.2,21);
%              >> [Z,R] = meshgrid(z,r);
%              >> I = [61.56  66.45  113.43  108.09 ];      
%              >> [Aphi] = B_Ellip_Cryogenic_170('aphi','cryogenic',I,R,Z);
%              >> figure(2); contour(Z,R,R.*Aphi);
%
% Magnetic field associated with GT170 GHz magnet, with currents Igun and Icav,
% as read on power supplies (i.e. a current -(Icav+Igun) is actually flowing in the 
% gun coil.
% The computation is based on the elliptic integral formulae that can be found in:
% Smythe: "Static and Dynamic Electricity", Third edition, McGraw-Hill, p.290.
%
%	Version :	1.0		J.-P. Hogge		Sep. 16, 2003
%               1.1     J.-P. Hogge     Jan. 12, 2005 Included 165Ghz FZK magnet, other options deactivated
%
%   mode    : 'bz','br','dbzbz','dbzdr','dbrdz','dbrdr', 'd2bzdz2' or 'aphi'
%   magnet	: 'fzk_165ghz_oi' only so far
%	     z	: a row vector (in [m]) or a 2D matrix containing the z locations where the field is to be evaluated
%	     r	: a row vector (in [m]) or a 2D matrix containing the r locations where the field is to be evaluated
%
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 140 GHz coils geometry [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        switch lower(magnet)
            
        case 'cryogenic', % Geometry of the 'as built' Cryogenic magnet for the EU 170GHz 1MW tube (Contract F4E OPE-552)
			 
			Ncoil   = 4;
            %                G1         G2         BC          MC
   		    Nturn   = [     959       2068       3898       26667      ];
 			rmin	= [0.181500   0.181450   0.147700    0.165300      ];
 			rmax	= [0.211360   0.212000   0.197200    0.240900      ];
 			zmin	= [0.084000   0.150400   0.207800    0.333110      ];
 			zmax    = [0.118650   0.185650   0.307850    0.637690      ];
			na		= 8. *[   10        10         10          10      ];   % Number of subcoils (radial direction
			nb		= 8. *[   10        10         10          10      ];   % Number of subcoils in longitudinal direction

		    Current = [   I(1)     -I(2)       -I(3)        I(4)       ];   % Current flowing in each power supply

  			rcen	= (rmin + rmax)/2; 					
  			zcen	= (zmin + zmax)/2;
  			dr		=  rmax - rmin; 
  			dz		=  zmax - zmin;
			coilsurf=  dr .* dz;
            
            JTot	=	Current .* Nturn;		% Total current flowing in each coil
		    J		=   JTot./(dr.*dz);     	% Current density in each coil
  	        Ic1     =   JTot ;                  % Total current flowing in each coil [A]

      case 'cryogenic_bestfit', % Geometry of the 'best fit' Cryogenic magnet for the EU 170GHz 1MW tube (Contract F4E OPE-552), 2018 05 25
			 
			Ncoil   = 4;
   		    Nturn   = [     959        2068       3898       26667      ];
 			rmin	= [0.186875    0.186540   0.146717    0.161725      ];
 			rmax	= [0.217632    0.202334   0.200795    0.242476      ];
 			zmin	= [0.0809181   0.153461   0.210388    0.332353      ];
 			zmax    = [0.115717    0.189897   0.304647    0.637898      ];
			na		= 8. *[   10        10         10          10      ];   % Number of subcoils (radial direction
			nb		= 8. *[   10        10         10          10      ];   % Number of subcoils in longitudinal direction

		    Current = [   I(1)     -I(2)       -I(3)        I(4)       ];   % Current flowing in each power supply

  			rcen	= (rmin + rmax)/2; 					
  			zcen	= (zmin + zmax)/2;
  			dr		=  rmax - rmin; 
  			dz		=  zmax - zmin;
			coilsurf=  dr .* dz;
            
            JTot	=	Current .* Nturn;		% Total current flowing in each coil
		    J		=   JTot./(dr.*dz);     	% Current density in each coil
  	        Ic1     =   JTot ;                  % Total current flowing in each coil [A]

        case 'design', % Geometry which appears in the SC CFT document (Contract EFDA/ 03-961)
			 
			Ncoil   = 7;
   		    Nturn   = [    3000       3000       4320        6200        9177        4566       7723  ];
 			rmin	= [0.135000   0.135000   0.135000    0.135000    0.190000    0.135000   0.170000  ];
 			rmax	= [0.145000   0.145000   0.200700    0.190000    0.237040    0.170000   0.204370  ];
 			zmin	= [0.078500   0.138500   0.183900    0.317000    0.317000    0.509000   0.509000  ];
 			zmax    = [0.098500   0.158500   0.282900    0.490000    0.490000    0.708000   0.708000  ];
			na		= 2. *[   10        10         10          10          10          10         10  ];   % Number of subcoils (radial direction
			nb		= 2. *[   10        10         10          10          10          10         10  ];   % Number of subcoils in longitudinal direction

		    Current = [   I(1)      I(2)        -I(3)        I(3)        I(3)        I(3)       I(3)  ];   % Current flowing in each power supply

  			rcen	= (rmin + rmax)/2; 					
  			zcen	= (zmin + zmax)/2;
  			dr		=  rmax - rmin; 
  			dz		=  zmax - zmin;
			coilsurf=  dr .* dz;
            
            JTot	=	Current .* Nturn;		% Total current flowing in each coil
		    J		=   JTot./(dr.*dz);     	% Current density in each coil
  	        Ic1     =   JTot ;                  % Total current flowing in each coil [A]
            
        case 'asg_modified_201206', 
%
% Geometry of the magnet as built with 2 layers removed on coils#3 and 3 on
% coil #5 after the problems with the conductor
% Option added on Decemebr 20th, 2006
%
% Reference current  
% [ 9.7 7.7 -96.2499 -96.2499 91.0124 91.0124 87.4723 87.4723 ] = [ I(1) I(2) -I(3) -I(3) I(4) I(4) I(5) I(5)] 
 
            Ncoil   = 8;
            Nturn   = [    218        218     2528.5       2626       7689    11517.5      6110    9656.5   ];
            rmin    = [0.13958    0.13958    0.13710    0.16952    0.13670    0.19101   0.13458   0.17197   ];
            rmax    = [0.142694   0.142694   0.16733    0.20005    0.18941    0.23666   0.16972   0.20412   ];
            zmin    = [0.07854    0.13854    0.18200    0.18200    0.31778    0.31748   0.50685   0.50685   ];
            zmax    = [0.09846    0.15846    0.28350    0.28350    0.48977    0.49017   0.71128   0.71128   ];
 
            rcen    = (rmin + rmax)/2;                  
            zcen    = (zmin + zmax)/2;
            dr      =  rmax - rmin; 
            dz      =  zmax - zmin;
            coilsurf=  dr .* dz;
 
            a       = rcen-dr/2; 
            b       = rcen+dr/2;
            l       = dz;
 
                        
            na      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils (radial direction
            nb      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils in longitudinal direction
 
            Current = [    I(1)       I(2)      -I(3)       -I(3)        I(4)        I(4)       I(5)       I(5)   ];   % Current flowing in each power supply
 
            JTot    =   Current .* Nturn;       % Total current flowing in each coil
            J       =   JTot./(dr.*dz);         % Current density in each coil
            Ic1     =   JTot ;                  % Total current flowing in each coil [A]
            
          case 'asg_report', 
%
% Geometry of the magnet as built with 2 layers removed on coils#3 and 3 on
% coil #5 after the problems with the conductor
% Option added on Decemebr 20th, 2006
%
% Reference current 
% [0  0  6.79541  1.967  88.2799]=[I(1) I(2) I(3) I(4) I(5)]

 
            Ncoil   = 8;
            Nturn   = [    218        218     2528.5       2626       7689    11517.5      6110    9656.5   ];
            rmin    = [0.13958    0.13958    0.13710    0.16952    0.13670    0.19101   0.13458   0.17197   ];
            rmax    = [0.142694   0.142694   0.16733    0.20005    0.18941    0.23666   0.16972   0.20412   ];
            zmin    = [0.07854    0.13854    0.18200    0.18200    0.31778    0.31748   0.50685   0.50685   ];
            zmax    = [0.09846    0.15846    0.28350    0.28350    0.48977    0.49017   0.71128   0.71128   ];
 
            rcen    = (rmin + rmax)/2;                  
            zcen    = (zmin + zmax)/2;
            dr      =  rmax - rmin; 
            dz      =  zmax - zmin;
            coilsurf=  dr .* dz;
 
            a       = rcen-dr/2; 
            b       = rcen+dr/2;
            l       = dz;
 
                        
            na      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils (radial direction
            nb      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils in longitudinal direction
 
            Current = [      I(1)       I(2)      -I(3)-I(5)  -I(3)-I(5)   I(4)+I(5)   I(4)+I(5)  I(5)       I(5)   ];   % Current flowing in each power supply
 
            JTot    =   Current .* Nturn;       % Total current flowing in each coil
            J       =   JTot./(dr.*dz);         % Current density in each coil
            Ic1     =   JTot ;                  % Total current flowing in each coil [A]
			
          otherwise, 
			disp('Unknown magnet')
	  			  
		end
	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct array of current loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Number of coils
	ncoil = length(rcen);

% Number of subcoils
	ntot  = sum(na.*nb);	

% Preallocate space
	r2    = zeros(1,ntot);
	z2    = zeros(1,ntot);
	Ic    = zeros(ntot,1);
	
	
% isub is the index of subcoil
% Positions of the coils subdivided
	isub = 0;	
	for icoil=1:ncoil
	   for ib=1:nb(icoil)
	      for ia=1:na(icoil)
				 isub 	  = isub +1;
				 r2(isub) = (rcen(icoil) - dr(icoil)/2) +  (dr(icoil)/nb(icoil))*(ib-0.5);
				 z2(isub) = (zcen(icoil) - dz(icoil)/2) +  (dz(icoil)/na(icoil))*(ia-0.5);
				 Ic(isub) = Ic1(icoil)            / (na(icoil)*nb(icoil));
           end
	   end
	end

%	printf('isub =',isub)


% Transform array of evaluation points in a one-dimensional array
%
	r1d = reshape(r,1,prod(size(r)));
	z1d = reshape(z,1,prod(size(z)));

%
% Suppress locations where Green function diverges (i.e. on the sources)
%
for i=1:length(r2),
  	   [I] = find((r1d.*r1d -r2(i)*r2(i) == 0) & (z1d.*z1d - z2(i)^2)==0);

	   if ~isempty(I)
	      r1d(I)=NaN;
	      z1d(I)=NaN;
	   end
	end
	
% Call greenem

	b = greenem_jph(mode,r1d,z1d,r2,z2);
    
    if(~iscell(mode))
        nbmodes=length({mode});
    else
        nbmodes=length(mode);
    end
    B=zeros(size(r,1),size(r,2),nbmodes);
	
% b is a rectangular matrix (length(r1d) x length(r2)) which needs to be multiplied by the 
% column vector Ic to get the proper contribution of each source.
	for i=1:size(b,3)
        temp = b(:,:,i) * Ic;
		
        B(:,:,i) = reshape(temp,size(r,1),size(r,2));
    end
	B(find((B==Inf)|(B==-Inf)))= 0;
	B(find(isnan(B)))	= 0;



	return
