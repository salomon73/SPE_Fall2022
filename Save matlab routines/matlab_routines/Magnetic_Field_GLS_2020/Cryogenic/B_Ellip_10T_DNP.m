function [B] = B_Ellip_10T_DNP(mode,magnet,I,r,z)
%
%    Call:     [B] = B_Ellip_10T_DNP(mode,magnet,r,z)
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
% 10T magnet for DNP gyrotron coils geometry [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(magnet)
    
    case 'new' % Geometry of the 10Tesla sup.cond. solenoid with 76 mm bore diameter
        
        Ncoil   = 13;
        Nturn   = [1670 2117    2960    4995    12793   59  57  56  54  59  57  56  54];
        rmin	= [52.592   60.069  76.919  82.502  90.129  103.638 104.187 104.735 105.283 103.638 104.187 104.735 105.283].*10^(-3);
        rmax	= [60.069   66.051  82.502  90.129  103.638 104.187 104.735 105.283 105.832 104.190  104.735 105.283 105.832].*10^(-3);
        zmin	= [-145 -145 -150 -150 -150 -150 -150 -150 -150 -150 -150 -150 -150].*10^(-3);
        zmax    = [145 145 150 150 150 150 150 150 150 150 150 150 150].*10^(-3);
        na		= 8. *[10 10 10 10 10 10 10 10 10 10 10 10 10];   % Number of subcoils (radial direction
        nb		= 8. *[10 10 10 10 10 10 10 10 10 10 10 10 10];   % Number of subcoils in longitudinal direction
        %   I = 31.31 A => B~ 1.7T at flange
        %   I = 32.35 A => B< 3T at z = 0 => 235 mm from top flange
        Current = [   I(1)     I(2)       I(3)        I(4) I(5) I(6) I(7) I(8) I(9) I(10) I(11) I(12) I(13)];   % Current flowing in each power supply
        
        rcen	= (rmin + rmax)/2;
        zcen	= (zmin + zmax)/2;
        dr		=  rmax - rmin;
        dz		=  zmax - zmin;
        coilsurf=  dr .* dz;
        
        JTot	=	Current .* Nturn;		% Total current flowing in each coil
        J		=   JTot./(dr.*dz);     	% Current density in each coil
        Ic1     =   JTot ;                  % Total current flowing in each coil [A]
        
    otherwise
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
