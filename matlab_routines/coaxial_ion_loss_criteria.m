%% wr+
H0=3.2e-14;
P0=8.66e-25;
r0=0.005;

qe=1.60217662E-19;
me=9.109383E-31;
eps_0=8.85418781762e-12;
vlight=299792458;

ne=6e17;
ne=5*logspace(10,19,5000);

kb=1.38e-23;
T=300;
P=1e-9*100000;
estimneutraldensity=P/(kb*T)

I=40;
V=90;
alpha=1.3;


B0=0.21;
L=0.48;
z=0.0;

% Davidson
prefix='Davidson';
phia=-60000;
phib=0;
R=1.5;
B0=0.21;
width=0.48;
b=0.16;
a=0.001;
r=0.0079;
zmin=0.0;
zmax=-width/2;
deltar=(R-1)/(R+1)*besseli(1,2*pi*r/width)*(width/2/pi)*(cos(2*pi*zmin/width)-cos(2*pi*zmax/width));
rb_=0.00764;
rbplus=0.0842;
Te=400; %eV
%1rbplus=0.020;

% Davidson test
prefix='Davidson';
phia=-60000;
phib=0;
R=1.5;
B0=0.21;
width=0.48;
b=0.16;
a=0.001;
rb_=0.0076;
r=9.7e-3 ;
zmin=0.0;
zmax=-width/8;
deltar=(R-1)/(R+1)*besseli(1,2*pi*r/width)*(width/2/pi)*(cos(2*pi*zmin/width)-cos(2*pi*zmax/width));
rbplus=rb_ + max(deltar,0.0026);
rcenter=r+deltar
Te=400; %eV
%1rbplus=0.020;


% rb_=0.01;
% rbplus=rb_+0.1;
% r=(rb_+rbplus)/2;
% r=0.015
% deltar=(R-1)/(R+1)*r;



%ITER
prefix='ITER' 
R=1.0001;
phia=-5000;
phib=0;
b=0.0792;
a=0.064;
r=0.080;
deltar=(R-1)/(R+1)*r;
rb_=0.0766;
rbplus=0.079;

spcharge=@(ne)-qe*ne/2/eps_0;

deltaphicloud=@(ne) Phi(r+deltar,rb_,rbplus,a,b,phia,phib, ne)-Phi(r,rb_,rbplus,a,b,phia,phib, ne);
deltaphivacuum=(phib-phia)/log(b/a)*log((r+deltar)/r);


figure
rforphi=linspace(a,b,1000);
phir=zeros(length(rforphi),1);
for i=1:length(rforphi)
    phir(i)=Phi(rforphi(i),rb_,rbplus,a,b,phia,phib, 6e17);
end
plot(rforphi,phir)
title('Phi')
figure
plot((rforphi(1:end-1)+rforphi(2:end))/2,-diff(phir)./diff(rforphi))
title('E_r')
Emaxcloud=1/(R-1)*deltaphicloud(ne);
Emaxvacuum=1/(R-1)*deltaphivacuum

rgrid=linspace(1e-3,16e-2,500);
zgrid=linspace(-0.24,0.24,1000);
GainedE=zeros(length(rgrid),length(zgrid));
for i=1:length(zgrid)
    drgrid=(R-1)/(R+1)*besseli(1,2*pi*rgrid/width)*(width/2/pi)*(cos(2*pi*0/width)-cos(2*pi*zgrid(i)/width));
    GainedE(:,i)=-(Phivacuum(rgrid,a,b,phia,phib)-Phivacuum(rgrid+drgrid,a,b,phia,phib));
end
figure
title(sprintf('\\phi_a=%3.2f kV \\phi_b=%3.2f kV R=%3.2f',phia/1e3,phib/1e3, R))
surface(zgrid,rgrid,GainedE);
xlabel('z')
ylabel('r')
c=colorbar;
c.Label.String='E [eV]';
meanE=mean(mean(abs(GainedE(20:375,ceil(length(zgrid)/2):floor(5*length(zgrid)/6))),2))/2
%%
figure('name',sprintf('%s Phi_a=%3.2f kV Phi_b=%3.2f kV r=%3.2f mm R=%3.2f',prefix,phia/1e3,phib/1e3,r*1e3, R));
loglog(ne,abs(Emaxcloud),'displayname','flat top density cloud')
hold on
loglog([ne(1) ne(end)], Emaxvacuum*[1 1],'displayname','Vacuum') 
xlabel('electron density [m^{-3}]')
ylabel('|E_{min}| [eV]')
title(sprintf('%s \\Phi_a=%3.2f kV \\Phi_b=%3.2f kV r=%3.2f mm R=%3.2f',prefix,phia/1e3,phib/1e3,r*1e3, R))
savefig(sprintf('%s_Phia%3.2f_Phib%3.2f_r%3.2f_R%3.2f.fig',prefix,phia/1e3,phib/1e3,r*1e3, R))


figure
subplot(2,2,1)
sgtitle(sprintf('%s \\Phi_a=%3.2f kV \\Phi_b=%3.2f kV r_0=%3.2f mm z_0=%3.2f mm R=%3.2f',prefix,phia/1e3,phib/1e3,(r+deltar)*1e3,z*1e3, R))
neinter=5e13;
vpar=linspace(1,vlight,1000);
vper=sqrt(-2*qe/me*deltaphicloud(neinter)/(R-1)+vpar.^2/(R-1))/vlight;
vpar=vpar/vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);

%Davidson vperp
phi0=Phi(r,rb_,rbplus,a,b,phia,phib, neinter);
vrz=sqrt(2/me*H0+2/me*qe*phi0-1/me^2*(P0/r+qe*Athet(r,z,B0,width,R))^2);

theta=2*pi*linspace(0,1,1000);
vthet=(P0/r+qe*Athet(r,z,B0,width,R))/me;
vperdav=sqrt((vrz^2*cos(theta).^2+vthet^2))/vlight;
vpardav=vrz*sin(theta)/vlight;

vpargauss=sqrt(Te*qe/me)*linspace(-1,1,200);
vpergauss=sqrt(2*Te*qe/me)*sqrt(1-vpargauss.^2*me/Te/qe);
vpargauss=[flip(vpargauss(2:end)) vpargauss]/vlight;
vpergauss=[-flip(vpergauss(2:end)) vpergauss]/vlight;

plot(vpar,vper,'b-')
hold on
plot(vpar,-vper,'b-')
plot(-vpar,vper,'b-')
plot(-vpar,-vper,'b-')
xlabel('\beta_z [m/s]')
ylabel('\beta_\perp [m/s]')
title(sprintf('electrons n_e=%3.2e m^{-3}',neinter))
grid
h(1)=plot(vpardav,vperdav,'r--','Displayname',sprintf('Dav H0=%#.2g eV P0=%#.2g kg m^2 s^{-1}',H0/qe,P0));
h(2)=plot(vpargauss,vpergauss,'g--','Displayname',sprintf('Gauss Te=%#.2g eV',Te));
legend(h)
legend('location','northoutside')
xlim([-1 1])
ylim([-1 1])

subplot(2,2,3)
neinter=5e17;
vpar=linspace(1,vlight,1000);
vper=sqrt(-2*qe/me*deltaphicloud(neinter)/(R-1)+vpar.^2/(R-1))/vlight;
vpar=vpar/vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);

%Davidson vperp
phi0=Phi(r,rb_,rbplus,a,b,phia,phib, neinter);
vrz=sqrt(2/me*H0+2/me*qe*phi0-1/me^2*(P0/r+qe*Athet(r,z,B0,width,R))^2);
theta=2*pi*linspace(0,1,1000);
vthet=(P0/r+qe*Athet(r,z,B0,width,R))/me;
vperdav=sqrt((vrz^2*cos(theta).^2+vthet^2))/vlight;
vpardav=vrz*sin(theta)/vlight;


vpargauss=sqrt(Te*qe/me)*linspace(-1,1,200);
vpergauss=sqrt(2*Te*qe/me)*sqrt(1-vpargauss.^2*me/Te/qe);
vpargauss=[flip(vpargauss(2:end)) vpargauss]/vlight;
vpergauss=[-flip(vpergauss(2:end)) vpergauss]/vlight;

plot(vpar,vper,'b-')
hold on
plot(vpar,-vper,'b-')
plot(-vpar,vper,'b-')
plot(-vpar,-vper,'b-')
xlabel('\beta_z [m/s]')
ylabel('\beta_\perp [m/s]')
title(sprintf('electrons n_e=%3.2e m^{-3}',neinter))
grid
h(1)=plot(vpardav,vperdav,'r--','Displayname',sprintf('Dav H0=%#.2g eV P0=%#.2g kg m^2 s^{-1}',H0/qe,P0));
h(2)=plot(vpargauss,vpergauss,'g--','Displayname',sprintf('Gauss Te=%#.2g eV',Te));
xlim([-1 1])
ylim([-1 1])

subplot(2,2,2)
neinter=5e13;
vpar=linspace(1,vlight,1000);
m=28.0134/1000*6.02214076e23;
vper=sqrt(2*qe/m*deltaphicloud(neinter)/(R-1)+vpar.^2/(R-1))/vlight;
vpar=vpar/vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
plot(vpar,vper,'b-')
hold on
plot(vpar,-vper,'b-')
plot(-vpar,vper,'b-')
plot(-vpar,-vper,'b-')
xlabel('\beta_z [m/s]')
ylabel('\beta_\perp [m/s]')
title(sprintf('N_2^+ ions n_e=%3.2e m^{-3}',neinter))
grid
xlim([-1 1])
ylim([-1 1])

subplot(2,2,4)
neinter=5e17;
vpar=linspace(1,vlight,1000);
vper=sqrt(2*qe/m*deltaphicloud(neinter)/(R-1)+vpar.^2/(R-1))/vlight;
vpar=vpar/vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
plot(vpar,vper,'b-')
hold on
plot(vpar,-vper,'b-')
plot(-vpar,vper,'b-')
plot(-vpar,-vper,'b-')
xlabel('\beta_z [m/s]')
ylabel('\beta_\perp [m/s]')
title(sprintf('N_2^+ ions n_e=%3.2e m^{-3}',neinter))
grid
savefig(sprintf('%s_Phia%3.2f_Phib%3.2f_r%3.2f_z%3.2f_R%3.2f_cone.fig',prefix,phia/1e3,phib/1e3,(r+deltar)*1e3,z*1e3, R))
xlim([-1 1])
ylim([-1 1])

%%
figure
r=double(M.rgrid);
r(r<1e-3)=1e-3;
r(r>16e-2)=16e-2;
phi=zeros(length(r),1);
for i=1:length(r)
    phi(i)=Phi(r(i),14e-2,15e-2,0.001,16e-2,-60000,0,0);
end
plot(r,phi)

phim=phi*ones(1,length(M.zgrid));

function Atheta=Athet(r,z,B0,width,R)
                Atheta=0.5*B0*(r-width/pi*(R-1)/(R+1)...
                .*besseli(1,2*pi*r/width).*cos(2*pi*z/width));
end

function phi=Phi(r,rbm,rbp,a,b,phia,phib, ne)
qe=1.60217662E-19;
me=9.109383E-31;
eps_0=8.85418781762e-12;
q=-qe;

phibp=1/log(b/a)*(ne/2*q/eps_0*rbm^2*log(b/rbp)*log(rbp/rbm)+(q*ne/2/eps_0*(rbp^2-rbm^2)*(log(rbp/a)-0.5)+phia)*log(b/rbp)+phib*log(rbp/a));
phibm=1/log(b/a)*(-q/eps_0*ne/2*rbm^2*log(rbm/a)*log(rbp/rbm)+(q/eps_0*ne/2*(rbp^2-rbm^2)*(log(b/rbp)+0.5)+phib)*log(rbm/a)+phia*log(b/rbm));
if(r<rbm && r>=a)
    phi=((phibm-phia)*log(r)+phia*log(rbm)-phibm*log(a))/log(rbm/a);
elseif(r>=rbm && r<=rbp)
    phi=-q/eps_0*ne/4*(r^2-rbp^2)+phibp+(phibp-phibm+q/eps_0*ne/4*(rbp^2-rbm^2))*log(r/rbp)/log(rbp/rbm);
elseif(r>rbp && r<= b)
    phi=((phib-phibp)*log(r)+phibp*log(b)-phib*log(rbp))/log(b/rbp);
else
    error('invalid radial position')
end

end
function phi=Phivacuum(r,a,b,phia,phib)
    phi=((phib-phia)*log(r)+phia*log(b)-phib*log(a))/log(b/a);
end