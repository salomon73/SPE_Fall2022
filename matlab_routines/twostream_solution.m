clear all
close all
%%

L=2*pi;             % spatial length of the system
DT=.1;              % time step
NT=1000;             % number of time steps
NG=128;              % number of grid points
N=2;             % number of particles
WP=1;               % plasma frequency (in normalized units)
QM=-1;              % charge of the particles (electrons)
V0=0.0;             % mean velocity
VT=0.0;             % thermal velocity
XP1=1;              % pertubation of positions
Q=WP^2/(QM*N/L);    % electron charge in normalized units
rho_back=-Q*N/L;    % background charge density (ions)
dx=L/NG;            % size of the grid cells

%% initial loading for the 2 Stream instability
xp=transpose(linspace(0,L-L/N,N));
vp=VT*randn(N,1);
pm=transpose([1:N]);pm=1-2*mod(pm,2);
vp=vp+pm.*V0;

% Perturbation
xp=xp+XP1*(L/N)*pm.*rand(N,1);
p=1:N;p=[p p];
un=ones(NG-1,1);
Poisson=spdiags([un -2*un un],[-1 0 1],NG-1,NG-1);

%% Prepare Diagnostics
time = 0:DT:(NT-1)*DT;
Phi=0; mat=0; Eg=0;
NFFT = 2^nextpow2(NG);
k = (2*pi/dx)/2*linspace(0,1,NFFT/2+1);
dv=V0/10; vlim=4*V0; vbin=-vlim:dv:vlim; NV=length(vbin);
mom = zeros(1,NT);
E_kin = zeros(1,NT);
E_pot = zeros(1,NT);
E_the = zeros(1,NT);
dens = zeros(NG,NT);
es_pot = zeros(NG,NT);
es_fld = zeros(NG,NT);
phi_fft = zeros(NG,NT);
vdf = zeros(NV,NT);

%% Main computational cycle
for it=1:NT
    
    %-----DIAGNOSTICS-----
    
    %total momentum
    
    mom(it) = sum(vp);
    E_kin(it) = 0.5*mean(vp.^2)*WP^2/QM^2;
    E_pot(it) = 0.5*sum(Eg.^2)/NG;
    E_the(it) = 0.5*mean((vp-mean(vp)).^2)*WP^2/QM^2;
    dens(:,it) = sum(mat)/dx;
    es_pot(:,it) = Phi;
    es_fld(:,it) = Eg;
    vdf(:,it) = histc(vp,vbin);
    phi_fft(:,it) = fft(Phi,NFFT)/NG;
    
    %---------------------
    
    % update xp and apply boundary conditions (bc)
    xp=mod(xp+vp*DT,L);
    
    % project particles to grid
    g1=floor(xp/dx);
    g=[g1;g1+1];
    % compute what fraction of the particle size that lies on the two nearest cells
    fraz1=1-abs(xp/dx-g1);
    fraz=[fraz1;1-fraz1];
    
    % apply bc on the projection
    out=(g<1);g(out)=g(out)+NG;
    out=(g>NG);g(out)=g(out)-NG;
    % use represention with a sparse matrix to save up memory
    mat=sparse(p,g,fraz,N,NG); 
    rho=transpose(full((Q/dx)*sum(mat))+rho_back); % charge density
    
    % computing fields
    Phi=Poisson\(-rho(1:NG-1)*dx^2);Phi=[Phi;0];
    Eg=([Phi(NG); Phi(1:NG-1)]-[Phi(2:NG);Phi(1)])/(2*dx);
    
    % projection q->p and update of vp
    vp=vp+mat*QM*Eg*DT;
    
    % plot phase space particle dynamics
    scatter(xp,vp,'.')
    %ylim([-V0*3 V0*3]);
    xlabel('X')
    ylabel('V')
    pause(0.01)
    
    
end

%%
figure;
plot(time,mom)
xlabel('time')
ylabel('momentum')

figure;
hold on
plot(time,E_kin,'k')
plot(time,E_pot,'r')
plot(time,E_the,'m')
plot(time,E_kin+E_pot,'b')
xlabel('time')
ylabel('Energy')

X=[0:dx:L-L/NG];
figure
surfc(time,X,dens(:,:)); shading interp; colorbar;
xlabel('time')
ylabel('X')
title('density')

figure
surfc(time,X,es_pot(:,:)); shading interp; colorbar;
xlabel('time')
ylabel('X')
title('ES potential')


figure
surfc(time,X,es_fld(:,:)); shading interp; colorbar;
xlabel('time')
ylabel('X')
title('ES field')


figure
surfc(time,k,abs(phi_fft(1:NFFT/2+1,:))); shading interp; colorbar;
xlabel('time')
ylabel('k')
title('ES potential (Fourier harmonics)')


figure
surfc(time,vbin,vdf); shading interp; colorbar;
xlabel('time')
ylabel('V')
title('f(v)')
