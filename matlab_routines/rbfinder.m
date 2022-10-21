%% wr+
H0=3.2e-14;
P0=8.66e-25;
r0=0.005;
consideredne=5e16;

%% wr-
% H0=2.0e-16;
% P0=-5.7e-25;
% r0=0.0025;

%% simulated
%          5e14     5e15     1e16     5e16     1e17     2e17     3e17     4e17      1e18    3e18     5e18     1e19
r_sim=   [0.00497, 0.0053,  0.0056,  0.0059,  0.00628, 0.0050,  0.0053,  0.0052, 0.0054,  0.0052,  0.0054];%,  0.0063];
rplussim=[0.0124,  0.0118,  0.0118,  0.0118,  0.0121,  0.0137,  0.0140,  0.016,   0.0163,  0.020,   0.025];%,   0.042];
nsim=    [4.8e14,  4.12e15, 7.5e15,  2.54e16, 4.19e16, 6.42e16, 8.55e16, 9.79e16, 1.34e17, 2.408e17, 2.416e17];%, 3.49e17];


B0=0.21;
Rcurv=1.5;
L=0.48;
z=0.0;

qe=1.60217662E-19;
me=9.109383E-31;
eps0=8.85418781762e-12;

ne=5*logspace(14,19,800);
rb_=zeros(size(ne));
rbplus=zeros(size(ne));
i=1;
minus=@(r) 1-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
%minus=@(r) 1+qe*((phib-phia)*log(r)+phia*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
format long;
rb_(i)=fzero(minus,r0);
plus=@(r) qe^2*ne(i)*rb_(i)/(eps0*H0)*(r-rb_(i)-rb_(i)*log(r/rb_(i)))+1-1/(2*me*H0)*(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2;
%plus=@(r) 1+qe*((phib-phia)*log(r)+phia*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
rbplus(i)=fzero(plus,5*rb_(i));
remainder=ones(size(ne));
remainder(i)=plus(rbplus(i));
for i=2:length(ne)
    try
        minus=@(r) 1-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
        %minus=@(r) 1+qe*((phib-phia)*log(r)+phia*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
        rb_(i)=fzero(minus,0.004);
        rb_t=rb_(i);
        plus=@(r) qe^2*ne(i)*rb_t/(eps0*H0)*(r-rb_t-rb_t*log(r/rb_t))+1-1/(2*me*H0)*(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2;
        %plus=@(r) 1+qe*((phib-phia)*log(r)+phia*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
        rbplus(i)=fzero(plus,1.5*rbplus(i-1));
        remainder(i)=plus(rbplus(i));
    catch
    end
end


f=figure;
%subplot(2,1,1)
semilogx(ne(remainder~=1),rbplus(remainder~=1)*1e3,'-','displayname','r_b^+')
hold on
semilogx(ne,rb_*1000,'-','displayname','r_b^-')
semilogx(nsim,rplussim*1e3,'+--','displayname','r_b^+ simulated','linewidth',1.5)
semilogx(nsim,r_sim*1e3,'+--','displayname','r_b^- simulated','linewidth',1.5)
legend
ylabel('r_b [mm]')
ylim([0 160])
xlabel('n_e(r_b^-) [m^{-3}]')
% subplot(2,1,2)
% semilogx(ne(remainder~=1),remainder(remainder~=1),'x--')
% ylabel('remainder [a.u.]')
% xlabel('n_e(r_b^-) [m^{-3}]')

f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
papsize=[10 10];
f.PaperSize=papsize;
print(f,'rbne','-dpdf','-fillpage')


i=find(consideredne>=ne,1,'last');
z=linspace(0,L/2,1000);
rb_z=zeros(size(z));
rbplusz=rb_z;
minus=@(r) 1-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z(1)/L)))^2/(2*me*H0);
rb_z(1)=fzero(minus,rb_(i));
rb_t=rb_z(1);
z_t=z(1);
plus=@(r) qe^2*ne(i)*rb_t/(eps0*H0)*(r-rb_t-rb_t*log(r/rb_t))+1-1/(2*me*H0)*(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z_t/L)))^2;
rbplusz(1)=fzero(plus,rbplus(i));
remainder(1)=plus(rbplusz(1));
for j=2:length(z)
minus=@(r) 1-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z(j)/L)))^2/(2*me*H0);
format long;
try
rb_z(j)=fzero(minus,rb_z(j-1));
rb_t=rb_z(j);
z_t=z(j);
plus=@(r) qe^2*ne(i)*rb_t/(eps0*H0)*(r-rb_t-rb_t*log(r/rb_t))+1-1/(2*me*H0)*(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z_t/L)))^2;
rbplusz(j)=fzero(plus,rbplusz(j-1));
remainder(j)=plus(rbplusz(j));
catch
end
end
figure()
z=[-flip(z(2:end)) z];
rb_z=[flip(rb_z(2:end)) rb_z];
rbplusz=[flip(rbplusz(2:end)) rbplusz];
plot(z(rb_z~=0),rb_z(rb_z~=0),'b')
hold on
plot(z(rbplusz~=0),rbplusz(rbplusz~=0),'r')
ylabel('r_b [m]')
ylim([0 0.16])
xlim([-L/2 L/2])
xlabel('z [m]')
title(sprintf('n_e=%7.2e [m^{-3}]',ne(i)))

ind=rbplusz~=0 & rb_z~=0& ~isnan(rbplusz) & ~isnan(rb_z);
ra=rb_z(ind);
rb=rbplusz(ind);
zvol=z(ind);
ra=(ra(1:end-1)+ra(2:end))/2;
rb=(rb(1:end-1)+rb(2:end))/2;
VOL=2*pi*diff(zvol).*min(ra).*(rb-ra);

nplasma=2116800
Nbesl=ne(i)*VOL;
Nbe=sum(Nbesl)
Nmacrosl=floor(Nbesl/Nbe*nplasma);
Nmacro=sum(Nmacrosl)
remain=nplasma-Nmacro
Nmacrosl(ceil(length(Nmacrosl)/2-remain/2):floor(length(Nmacrosl)/2+remain/2))=Nmacrosl(ceil(length(Nmacrosl)/2-remain/2):floor(length(Nmacrosl)/2+remain/2))+1;
Nmacrosl(floor(length(Nmacrosl)/2))=Nmacrosl(floor(length(Nmacrosl)/2))+nplasma-sum(Nmacrosl);
Nmacro=sum(Nmacrosl)
remain=nplasma-Nmacro


msim=me*Nbe/nplasma
qsim=qe*Nbe/nplasma

datas.ra=ra;
datas.rb=rb;
datas.z=zvol;
datas.npartsslice=Nmacrosl;
datas.m=me;
datas.weight=sum(Nbesl)/nplasma;
datas.q=-qe;
datas.H0=H0;
datas.P0=P0;
datas.temperature=10000;
datas.velocitytype=2;
datas.radialtype=1;