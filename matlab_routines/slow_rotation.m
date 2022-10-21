B0=0.21;
me=9.109e-31;
qe=1.6e-19;
omegace=qe*B0/me;
R=1.5;
z=0.0;
L=0.48;
invr=1-(R-1)/(R+1)*cos(2*pi*z/L);

H0=3.2e-14;
P0=-1.66e-23;

H0=1.2e-14
P0=-5.66e-25


r0=sqrt(4*H0/me/omegace^2);
rplus=r0/invr*sqrt(1-P0*omegace/2/H0*invr+sqrt(1-P0*omegace/H0*invr))
rminus=r0/invr*sqrt(1-P0*omegace/2/H0*invr-sqrt(1-P0*omegace/H0*invr))


r=(rplus+rminus)/2
Wtheta=(P0/me/r+qe*invr/2*B0*r/me)/r
1/2*me*(Wtheta*r)^2
omegace

zmax=L/2/pi*acos((R+1)/(R-1)*(1-H0/P0/omegace))

H0=P0*omegace*invr
