%% wr+
H0=3.2e-14;
P0=8.66e-25;
r0=0.005;
ne=5e14;

B0=0.21;
Rcurv=1.5;
L=0.48;
z=0.0;
consideredphia=-60000;

qe=1.60217662E-19;
me=9.109383E-31;
eps0=8.85418781762e-12;
phia=-linspace(0,90000,80);
phib=0;
b=0.16;
a=0.001;

%% simulated

r_sim=   [5.8e-3,  5.4e-3];%,  0.0063];
rplussim=[1.11e-2, 1.20e-2];%,   0.042];
phiasim= [-45000,  -30000];%, 3.49e17];



rb_phi=zeros(size(phia));
rbplusphi=zeros(size(phia));
i=1;
minus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
format long;
rb_phi(i)=fzero(minus,r0);
%plus=@(r) qe^2*ne(i)*rb_(i)/(eps0*H0)*(r-rb_(i)-rb_(i)*log(r/rb_(i)))+1-1/(2*me*H0)*(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2;
plus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
rbplusphi(i)=fzero(plus,5*rb_phi(i));
remainderphi=ones(size(phia));
remainderphi(i)=plus(rbplusphi(i));
for i=2:length(phia)
    try
        minus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
        rb_phi(i)=fzero(minus,rb_phi(i-1));
        rb_t=rb_phi(i);
        plus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z/L)))^2/(2*me*H0);
        rbplusphi(i)=fzero(plus,rbplusphi(i-1));
        remainderphi(i)=plus(rbplusphi(i));
    catch
    end
end

f=figure;
%subplot(2,1,1)
plot(phia(remainderphi~=1)/1e3,rb_phi(remainderphi~=1)*1e3,'x-','displayname','r_b^-')
hold on
plot(phia(remainderphi~=1)/1e3,rbplusphi(remainderphi~=1)*1e3,'x-','displayname','r_b^+')
plot(phiasim/1e3,rplussim*1e3,'+--','displayname','r_b^+ simulated','linewidth',1.5)
plot(phiasim/1e3,r_sim*1e3,'+--','displayname','r_b^- simulated','linewidth',1.5)
legend
ylabel('r_b [mm]')
ylim([0 160])
xlabel('\Phi_a [kV]')
% subplot(2,1,2)
% semilogx(ne(remainder~=1),remainder(remainder~=1),'x--')
% ylabel('remainder [a.u.]')
% xlabel('n_e(r_b^-) [m^{-3}]')

f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
papsize=[10 10];
f.PaperSize=papsize;
print(f,'rbne','-dpdf','-fillpage')



i=find(consideredphia>=phia,1,'first');
z=linspace(0,L/2,6000);
rb_z=zeros(size(z));
rbplusz=rb_z;

minus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z(1)/L)))^2/(2*me*H0);
rb_z(1)=fzero(minus,rb_phi(i));
rb_t=rb_z(1);
z_t=z(1);
plus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z(1)/L)))^2/(2*me*H0);        
rbplusz(1)=fzero(plus,rbplusphi(i));
remainder(1)=plus(rbplusz(1));
for j=2:length(z)
minus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z(j)/L)))^2/(2*me*H0);
format long;
try
rb_z(j)=fzero(minus,rb_z(j-1));
z_t=z(j);
plus=@(r) 1+qe*((phib-phia(i))*log(r)+phia(i)*log(b)-phib*log(a))/(log(b/a))/H0-(P0/r+qe*0.5*B0*(r-L/pi*(Rcurv-1)/(Rcurv+1)*besseli(1,2*pi*r/L)*cos(2*pi*z_t/L)))^2/(2*me*H0);        
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
plot(z(rb_z~=0),rbplusz(rb_z~=0),'r')
ylabel('r_b [m]')
ylim([0 0.16])
xlim([-L/2 L/2])
xlabel('z [m]')
title(sprintf('\\phi_a=%5.2f kV',phia(i)/1e3))

ind=rbplusz~=0 & rb_z~=0& ~isnan(rbplusz) & ~isnan(rb_z);
ra=rb_z(ind);
rb=rbplusz(ind);
zvol=z(ind);
ra=(ra(1:end-1)+ra(2:end))/2;
rb=(rb(1:end-1)+rb(2:end))/2;
VOL=2*pi*diff(zvol).*min(ra).*(rb-ra);

nplasma=2116800
Nbesl=ne*VOL;
Nbe=sum(Nbesl)
Nmacrosl=floor(Nbesl/Nbe*nplasma);
Nmacro=sum(Nmacrosl)
remain=nplasma-Nmacro
Nmacrosl(ceil(length(Nmacrosl)/2-remain/2):floor(length(Nmacrosl)/2+remain/2))=Nmacrosl(ceil(length(Nmacrosl)/2-remain/2):floor(length(Nmacrosl)/2+remain/2))+1;
Nmacrosl(floor(length(Nmacrosl)/2))=Nmacrosl(floor(length(Nmacrosl)/2))+nplasma-sum(Nmacrosl);
Nmacro=sum(Nmacrosl)
remain=nplasma-Nmacro


msim=me*Nbe/nplasma
qsim=-qe*Nbe/nplasma

datas.ra=ra;
datas.rb=rb;
datas.z=zvol;
datas.npartsslice=Nmacrosl;
datas.m=me;
datas.weight=Nbe/nplasma;
datas.q=-qe;
datas.H0=H0;
datas.P0=P0;
datas.temperature=10000;
datas.velocitytype=2;
datas.radialtype=1;