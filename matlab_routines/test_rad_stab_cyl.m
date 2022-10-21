a=1e-3;
b=0.03;
phia=-5e3;
phib=0;
rm=2.8e-3;
rp=0.02;%rm+3e-3;

nB=200;
nne=500;

B=linspace(0.01,3,nB);
ne=logspace(14,20,nne);

vlight=299792458;
qe=1.60217662E-19;
me=9.109383E-31;
eps_0=8.85418781762E-12;
kb=1.38064852E-23;

q=-qe;
C=-q*ne/(4*eps_0)*(rp^2-rm^2)/log(rp/rm)...
    +(phib-phia)/log(b/a)...
    +(q*ne/(2*eps_0)*rm^2*log(b*rm/(a*rp))...
   +q*ne/(2*eps_0)*(rp^2-rm^2)*(log(b/rp)+0.5*(1+log(a/b)/log(rp/rm))))/log(b/a);
C=1/2*q*ne*rm^2/eps_0;
omegape2=qe^2*ne/eps_0/me;
omegace2=qe^2*B.^2/me^2;

r=linspace(rm,rp,100);
for i=1:nB
    for j=1:nne
        [stab1(i,j),id]=min(1-2*omegape2(j)/omegace2(i)+4*q*C(j)/me/omegace2(i)./r.^2);
        rmin(i,j)=r(id);
    end
end

figure
contourf(B,ne,stab1'>0,'edgecolor','none')
set(gca,'YScale','log')
title(sprintf('r_m=%2.1fmm r_p=%2.1fmm, \\phi_a=%2.1fkV, \\phi_b=%2.1fkV',rm*1e3,rp*1e3,phia/1e3,phib/1e3))

