%%
rAthetpos=30;

fieldstep=length(M.t2d);%-100:5:length(M.t2d);
Vpar=M.fluidUR(:,:,fieldstep).*(M.Br./M.B)' + (M.Bz./M.B)'.*M.fluidUZ(:,:,fieldstep);
N=M.N(:,:,fieldstep);
Gammapar=mean(N.*Vpar,3);


levels=M.rAthet(rAthetpos,floor(M.nz/2)+1)*[1 1];

rAthetposend=find(levels(1)>=M.rAthet(:,1),1,'last');
rleft=sqrt(levels(1)/(M.rAthet(rAthetposend,1)/M.rgrid(rAthetposend)^2));

thefig=figure;
sp(1)=subplot(1,3,1);
contourf(M.zgrid,M.rgrid,mean(N,3),'edgecolor','none');
hold on;
contour(M.zgrid,M.rgrid,M.rAthet,levels,'r-','linewidth',1.5,'Displayname','Magnetic field lines');
c=colorbar;
c.Label.String='n_e [m^{-3}]';
sp(2)=subplot(1,3,2);
contourf(M.zgrid,M.rgrid,Gammapar,'edgecolor','none');
hold on
contour(M.zgrid,M.rgrid,M.rAthet,levels,'r-','linewidth',1.5,'Displayname','Magnetic field lines');
c=colorbar;
c.Label.String='\Gamma_{par} [m^{-2}s^{-1}]';
sp(3)=subplot(1,3,3);
contourf(M.zgrid,M.rgrid,mean(Vpar,3),'edgecolor','none');
hold on
contour(M.zgrid,M.rgrid,M.rAthet,levels,'r-','linewidth',1.5,'Displayname','Magnetic field line');
c=colorbar;
c.Label.String='V_{par} [ms^{-1}]';
linkaxes(sp);
xlabel(sp,'z [m]')
ylabel(sp,'r [m]')
ylim(sp,[M.rgrid(1) M.rgrid(sum(M.nnr(1:2)))])
M.savegraph(thefig,sprintf('%s/%s_finalwellfluxline',M.folder,M.name),[12,10])


%%
dN=1e13*M.weight;
dN=0;
rp=9.85e-3;
rm=7e-3;
V=2*pi*(rp^2-rm^2);
n2=mean(mean(mean(N(rAthetposend+(0:1),[1 end],:),3)));
vpar2=(dN/V/n2)^2;
actvpar2=mean(mean(abs(mean(Vpar(rAthetposend+(0:1),[1 end],:),3))))^2;
b=M.rgrid(end);
a=M.rgrid(1);
vd1=((M.potout-M.potinn)*M.phinorm/M.rgrid(rAthetpos))/M.Bz(floor(M.nz/2)+1,rAthetpos)/log(b/a);
%vd1=-M.Er(rAthetpos+5,M.nz/2+1,end)/M.Bz(M.nz/2+1,rAthetpos+5)
vinit=sqrt(M.kb*10000/M.me)
%vinit=sqrt(80*M.qe/M.me)
vd2=((M.potout-M.potinn)*M.phinorm/rleft)/M.Bz(1,rAthetposend)/log(b/a);
deltaphi=-0.5*M.me/M.qe*(vpar2-(vd1+vinit)^2*(1-M.Rcurv))
ratioparper=vpar2/((vd1+vinit)^2*(1-M.Rcurv))

%%


thefig=figure;
%fieldstep=fieldstep(end);
plotaxes(1)=subplot(1,2,1,'Parent',thefig);
plotaxes(2)=subplot(1,2,2,'Parent',thefig);
dens=mean(M.N(:,:,fieldstep),3);
model=M.potentialwellmodel(fieldstep);
z=model.z;
r=model.r;
pot=model.pot;
rathet=model.rathet;
[Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,floor(M.nz/2)+1));
densplot=zeros(size(Zmesh,1),size(Zmesh,2),length(fieldstep));
for i=1:length(fieldstep)
densplot(:,:,i)=griddata(Zmesh,M.rAthet,dens,Zmesh,Rmesh,'natural');
end
plot(plotaxes(1),M.zgrid,mean(densplot(rAthetpos,1:end,:),3));
xlim(plotaxes(1),[M.zgrid(1) M.zgrid(end)])
xlabel(plotaxes(1),'z [m]')
ylabel(plotaxes(1),'n [m^{-3}]')
title(plotaxes(1),'Density')
plotpot=zeros(size(Zmesh,1),size(Zmesh,2),length(fieldstep));
for i=1:length(fieldstep)
plotpot(:,:,i)=griddata(z,rathet,pot(:,i),Zmesh,Rmesh,'natural');
end
plotpot=mean(plotpot,3);
plot(plotaxes(2),M.zgrid(1:end),plotpot(rAthetpos,1:end));
hold on
plot(plotaxes(2),M.zgrid(1:end),deltaphi*ones(size(M.zgrid)),'k--','displayname','Analytical')
xlim(plotaxes(2),[M.zgrid(1) M.zgrid(end)])
xlabel(plotaxes(2),'z [m]')
ylabel(plotaxes(2),'depth [eV]')
title(plotaxes(2),'Well')
sgtitle(thefig,sprintf('r(0)=%1.2f [mm] t= %1.2f [ns]',M.rgrid(rAthetpos)*1e3,M.t2d(fieldstep(end))*1e9))

M.savegraph(thefig,sprintf('%s/%s_finalwell',M.folder,M.name),[15,10])




