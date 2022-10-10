%%
Pn=1e-4;%mbar
kb=1.38064852e-23;
T=300;
estimneutraldensity=Pn*100/(kb*T);
nb=estimneutraldensity;

qe=1.60217662E-19;
me=9.109383E-31;
eps_0=8.85418781762e-12;
vlight=299792458;

%%
Ee=[500 1000]; %eV
gridsize1=ceil(sqrt(length(Ee)));
gridsize2=ceil(length(Ee)/gridsize1);

fighandle=figure;
for i=1:length(Ee)
subplot(gridsize2,gridsize1,i)
ve=sqrt(Ee(i)*qe/me); %m/s
gamma=1/(1-ve^2/vlight^2);
ne=5*logspace(13,17,3000);
B0=0.28;
sigion=sigmabeb(15.58,54.91,2,1,Ee(i));
sigion=sigmabeb(41.72,71.13,2,1,Ee(i))+sigion;
sigion=sigmabeb(21,63.18,2,1,Ee(i))+sigion;
sigion=sigmabeb(17.07,44.30,4,1,Ee(i))+sigion;
sigion=M.sigio(Ee(i));

weight=3.83e6;
ncreated=5e17*nb*sigion*ve*(0.030*2*pi*(0.078^2-0.076^2))
nbsimucr=ncreated/weight
Coulomblog=log(sqrt(eps_0*Ee(i)/qe*ne)/qe^2*2*pi*eps_0*me*ve^2);

ne_0=1;
tauion=1./(nb*sigion*ve);
tauloss=tauion*log(ne/ne_0);
taucycl=2*pi*me/qe/B0;
taupe=2*pi*sqrt(eps_0*me/qe^2./ne);
tauthermal=1./(ne*qe^4.*Coulomblog./(me^2*ve^3*4*pi*eps_0^2));

ldw=1.5;
loglog(ne,tauion*ones(size(ne)),'displayname','\tau_{ion}','linewidth',ldw)
hold on
loglog(ne,tauloss,'displayname','t_{accum}','linewidth',ldw)
loglog(ne,tauthermal,'displayname','\tau_{thermal}','linewidth',ldw)
loglog(ne,taucycl*ones(size(ne)),'displayname','\tau_{ce}','linewidth',ldw)
loglog(ne,taupe,'displayname','\tau_{pe}','linewidth',ldw)
legend('location','eastoutside')
%title(sprintf('E_e=%#.3g eV',Te(i)))
xlabel('n_e [m^{-3}]')
ylabel('\tau [s]')
grid on
end
sgtitle(sprintf('E_e=%3.2g [eV], P_n=%3.2g [mbar]',Ee(1),Pn))

%% Saves the given figure as a pdf and a .fig 
fighandle.PaperUnits='centimeters';
papsize=[12 10];
fighandle.PaperSize=papsize;
name=sprintf('timescaleE%2.0eP%2.0embar',Ee(1),Pn);
print(fighandle,name,'-dpdf','-fillpage')
savefig(fighandle,name)

% Te=500;
% ve=sqrt(Te*qe/me); %m/s
% gamma=1/(1-ve^2/vlight^2);
% time=logspace(-5,1,500);
% dens=exp(time/tauion);
% Coulomblog=log(sqrt(eps_0*Te/qe*dens)/qe^2*2*pi*eps_0*me*ve^2);
% tauthermal=1./(ve*dens*qe^4.*Coulomblog./(gamma^2*me^2*ve^4*2*pi*eps_0^2));
% figure
% yyaxis left
% loglog(time,dens)
% yyaxis right
% loglog(time,tauthermal)




% Sigma BEB see (https://physics.nist.gov/PhysRefData/Ionization/intro.html)
% B binding energy 
% U average kinetic energy
% N electron occupation number
% Q dipole constant
% T incident energy
function sig=sigmabeb(B,U,N,Q,T)
t=T/B; %incident energy
u=U/B; % Average kinetic energy
% B:
a0=0.52918e-10;
R=13.6057;
S=4*pi*a0^2*N*(R/B)^2;
n=1;
sig=S/(t+(u+1)/n)*(Q*log(t)*0.5*(1-t^-2)+(2-Q)*(1-t^-1-log(t)/(t+1)));
end

