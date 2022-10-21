%    Computes the magnetic field of each submagnet in the asg
%    superconducting magnet assembly
%
%   The magnetic field is calculated for a unit current at positions
%   defined by r and z.
%   The fields are then saved in asg_magnet.mat

%   Evaluation grid space

z=linspace(-0.3,0.3,6000);
r=linspace(0.0,0.1,1000);
% z=linspace(-0.15,0.3,100);
% r=linspace(0.0,0.1,100);
[Z,R] = meshgrid(z,r);

mode={'aphi','bz','br'};


% Transform array of evaluation points in a one-dimensional array
%
r1d = reshape(R,1,numel(R));
z1d = reshape(Z,1,numel(Z));

%
% Geometry of the asg magnet as built with 2 layers removed on coils#3 and 3 on
% coil #5 after the problems with the conductor
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

% Current = [      I(1)       I(2)      -I(3)-I(5)  -I(3)-I(5)   I(4)+I(5)   I(4)+I(5)  I(5)       I(5)   ];   % Current flowing in each power supply
%
JTot    =   Nturn;       % Total current flowing in each coil
% J       =   JTot./(dr.*dz);         % Current density in each coil
Ic1     =   JTot ;                  % Total current flowing in each coil [A]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct array of current loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Number of coils
ncoil = length(rcen);

if(~iscell(mode))
    nbmodes=length({mode});
else
    nbmodes=length(mode);
end
subcoils=cell(ncoil,nbmodes);
for icoil=1:ncoil
    for imode=1:nbmodes
        subcoils{icoil,imode}=zeros(length(r),length(z));
    end
end
for icoil=1:ncoil
    %subcoils{icoil}=zeros(length(r),length(z),nbmodes);
    B=zeros(length(r),length(z),nbmodes);
    parfor ir=1:length(r)
        [Z,R] = meshgrid(z,r(ir));
        r1d = reshape(R,1,numel(R));
        z1d = reshape(Z,1,numel(Z));
        
        % Number of subcoils
        ntot  = sum(na(icoil).*nb(icoil));
        
        % Preallocate space
        r2    = zeros(1,ntot);
        z2    = zeros(1,ntot);
        Ic    = zeros(ntot,1);
        
        
        % isub is the index of subcoil
        % Positions of the coils subdivided
        isub = 0;
        for ib=1:nb(icoil)
            for ia=1:na(icoil)
                isub 	  = isub +1;
                r2(isub) = (rcen(icoil) - dr(icoil)/2) +  (dr(icoil)/nb(icoil))*(ib-0.5);
                z2(isub) = (zcen(icoil) - dz(icoil)/2) +  (dz(icoil)/na(icoil))*(ia-0.5);
                Ic(isub) = Ic1(icoil)            / (na(icoil)*nb(icoil));
            end
        end
        
        
        
        
        %
        % Suppress locations where Green function diverges (i.e. on the sources)
        %
        for i=1:length(r2)
            [I] = find((r1d.*r1d -r2(i)*r2(i) == 0) & (z1d.*z1d - z2(i)^2)==0);
            
            if ~isempty(I)
                r1d(I)=NaN;
                z1d(I)=NaN;
            end
        end
        
        % Call greenem
        
        b = greenem_jph(mode,r1d,z1d,r2,z2);
        
        
        temp2=zeros(size(R,2),nbmodes);
        
        % b is a rectangular matrix (length(r1d) x length(r2)) which needs to be multiplied by the
        % column vector Ic to get the proper contribution of each source.
        for i=1:size(b,3)
            temp = b(:,:,i)*Ic;
            temp((temp==Inf)|(temp==-Inf))= 0;
            temp(isnan(temp))	= 0;
            temp2(:,i) = reshape(temp,size(R,1),size(R,2));
        end
        B(ir,:,:)=temp2(:,:);
        
    end
    for imode=1:nbmodes
        subcoils{icoil,imode}=B(:,:,imode);
    end
end
save('asg_magnet_red.mat','r','z','subcoils','mode','Ncoil','na','nb','Nturn','rmin','rmax','zmin','zmax','-v7.3')


