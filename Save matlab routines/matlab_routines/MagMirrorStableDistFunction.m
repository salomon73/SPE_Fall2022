% Physical constants
vlight=299792458;
qe=1.60217662E-19;
me=9.109383E-31;
eps_0=8.85418781762E-12;
kb=1.38064852E-23;

H0=200000*1.6e-19;
B0=0.21;

plasmadim=[0, 0.16, 0.00749, 0.0075];
n0=-5e14;
nplasma=100000
%n0=-4e16

qsim=pi*(plasmadim(2)-plasmadim(1))*(plasmadim(4)^2-plasmadim(3)^2)*n0*qe/nplasma 
msim=abs(qsim)/qe*me

omegace=qe*B0/me;
omegape=sqrt(qe*qe*abs(n0)/eps_0/me);

r0=sqrt(2*H0/me/omegace^2)

vperp=omegace*r0/vlight


H0=200000*1.6e-19
Rcurv=1.5
width=0.64

deltar=sqrt((Rcurv-1)/(Rcurv+1))
rb=r0*(1+deltar)/(1-(Rcurv-1)/(Rcurv+1))
ra=r0*(1-deltar)/(1-(Rcurv-1)/(Rcurv+1))

v=[0.5*omegace*r0^2/ra,0.5*omegace*r0^2/rb]/vlight


%%
% Computation of the axial trapped resonance frequency

r=0.05:0.001:0.06;
z=-0.035:0.001:0.035;
B0=0.2;
Er=-30000;
L=0.05;
% Br=-B0*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*sin(2*pi*z/L)
% Bz=B0*(1-(Rcurv-1)/(Rcurv+1)*besseli(0,2*pi*r/L)*cos(2*pi*z/L))
% B=sqrt(Br.^2+Bz.^2)
