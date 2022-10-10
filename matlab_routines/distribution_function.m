%% Show the particles velocity histogram at the position set using ginput
% The histogram is compared between timesteppart and end

timesteppart=length(M.tpart);

fig=figure;
ax1=gca;

timestepN=find(M.tpart(end)==M.t2d,1);
Ndistrib=M.N(:,:,timestepN);
h=contourf(ax1,M.zgrid,M.rgrid,Ndistrib);
%set(h, 'EdgeColor', 'none');
hold on
[r,z]=find(Ndistrib~=0);
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= 'n [m^{-3}]';
view(ax1,2)
%set(ax1,'colorscale','log')


% [x,y]=ginput(1);
% Zindex=find(x>M.zgrid,1,'last');
% Rindex=find(y>M.rgrid,1,'last');
%  Rindex=155;
%  Zindex=344;

Z=(M.zgrid(Zindex));
R=(M.rgrid(Rindex));

plot(Z,R,'rx','Markersize',12);


% Rindex=16; Zindex=64;
% Rindex=28; Zindex=floor(size(M.zgrid,1)/2)+1;
nbp=min(M.R.nparts,M.nbparts(1));
Rp=M.R(1:nbp,1,false);
Zp=M.Z(1:nbp,1,false);
Indices=Rp>=M.rgrid(Rindex) & Rp<M.rgrid(Rindex+1) & Zp>=M.zgrid(Zindex) & Zp<M.zgrid(Zindex+1);

Vr=M.VR(Indices,1,false)/M.vlight;
Vz=M.VZ(Indices,1,false)/M.vlight;
Vthet=M.VTHET(Indices,1,false)/M.vlight;
Vperp=sqrt(Vthet.^2+Vr.^2);
Rp=M.R(Indices,1,false);
Zp=M.Z(Indices,1,false);
Thetp=M.THET(Indices,1,false);

nbp=min(M.R.nparts,M.nbparts(timesteppart));
Rend=M.R(1:nbp,timesteppart,false);
Zend=M.Z(1:nbp,timesteppart,false);
Indicesend=Rend>=M.rgrid(Rindex) & Rend<M.rgrid(Rindex+1) & Zend>=M.zgrid(Zindex) & Zend<M.zgrid(Zindex+1);

Vrend=M.VR(Indicesend,timesteppart,false)/M.vlight;
Vzend=M.VZ(Indicesend,timesteppart,false)/M.vlight;
Vthetend=M.VTHET(Indicesend,timesteppart,false)/M.vlight;
Vperpend=sqrt(Vthetend.^2+Vrend.^2);
Rend=M.R(Indicesend,timesteppart,false);
Zend=M.Z(Indicesend,timesteppart,false);
Thetend=M.THET(Indicesend,timesteppart,false);

binwidth=abs(max(Vrend)-min(Vrend))/sqrt(length(Vrend));
f=figure('Name',sprintf("%s v distrib",M.file));
tstudied=0;
legtext=sprintf("t=%2.1f [ns]",M.tpart(end)*1e9);
ax1=subplot(1,3,1);
if(length(Vr)>1)
h1=histogram(Vr,'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
%h1=histfit(ax1,Vr);
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
end
hold on
h1=histogram(Vrend,'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
%h1=histfit(ax1,Vrend);
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
ylabel('counts')
xlabel('\beta_r')

grid on

ax2=subplot(1,3,2);
binwidth=abs(max(Vthetend)-min(Vthetend))/sqrt(length(Vthetend));
if(length(Vthet)>1)
h1=histogram(Vthet(:,1),'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
end
%h1=histfit(ax2,Vthet(:,1));
hold on
h1=histogram(Vthetend,'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
%h1=histfit(ax2,Vthetend);
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
ylabel('counts')
xlabel('\beta_\theta')
legend(ax2,'location','northoutside','orientation','vertical')
grid on

ax3=subplot(1,3,3);
%h1=histogram(Vz(:,1),'Binwidth',binwidth,'Normalization','probability','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
binwidth=abs(max(Vzend)-min(Vzend))/sqrt(length(Vzend));
if(length(Vz)>1)
h1=histogram(ax3,Vz(:,1),'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
end
hold on
%h1=histogram(Vzend,'Binwidth',binwidth,'Normalization','probability','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
h1=histogram(ax3,Vzend,'Binwidth',binwidth,'Normalization','count','DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
set(h1,'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(end)*1e9));
ylabel('counts')
xlabel('\beta_z')
grid on
f=gcf;
sgtitle(sprintf('R=%1.2e[m] Z=%1.2e[m] dt=%1.2e[ns]',M.rgrid(Rindex+1),M.zgrid(Zindex+1),M.dt*1e9))
f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
papsize=[16 10];
f.PaperSize=papsize;
[~, name, ~] = fileparts(M.file);
print(f,sprintf('%sParts_V_RZ',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sParts_V_RZ',name))


%% Shows the particles phase space at position (rindex, zindex) 
% and times tinit and timesteppart
timeinit=1;
timestepNinit=find(M.tpart(timeinit)==M.t2d);
timesteppart=length(M.tpart);
timestepNend=find(M.tpart(timesteppart)==M.t2d);

Rp=M.R(:,timeinit,false);
Zp=M.Z(:,timeinit,false);
Indices=Rp>=M.rgrid(Rindex) & Rp<M.rgrid(Rindex+1) & Zp>=M.zgrid(Zindex) & Zp<M.zgrid(Zindex+1);

costhet=M.Bz(Zindex,Rindex)/M.B(Zindex,Rindex);
sinthet=M.Br(Zindex,Rindex)/M.B(Zindex,Rindex);

Vr=M.VR(Indices,timeinit,false)/M.vlight;
Vz=M.VZ(Indices,timeinit,false)/M.vlight;
Vthet=M.VTHET(Indices,timeinit,false)/M.vlight;
Vperp=M.Vperp(Indices,timeinit,false)/M.vlight;
Vpar=M.Vpar(Indices,timeinit,false)/M.vlight;
Rp=M.R(Indices,timeinit,false);
Zp=M.Z(Indices,timeinit,false);
Thetp=M.THET(Indices,timeinit,false);

Rend=M.R(:,timesteppart,false);
Zend=M.Z(:,timesteppart,false);
Indicesend=Rend>=M.rgrid(Rindex) & Rend<M.rgrid(Rindex+1) & Zend>=M.zgrid(Zindex) & Zend<M.zgrid(Zindex+1);

Vrend=M.VR(Indicesend,timesteppart,false)/M.vlight;
Vzend=M.VZ(Indicesend,timesteppart,false)/M.vlight;
Vthetend=M.VTHET(Indicesend,timesteppart,false)/M.vlight;
Vperpend=M.Vperp(Indicesend,timesteppart,false)/M.vlight;
Vparend=M.Vpar(Indicesend,timesteppart,false)/M.vlight;
Rend=Rend(Indicesend);
Zend=Zend(Indicesend);
Thetend=M.THET(Indicesend,timesteppart,false);

legendinit=sprintf('t=%1.3g [s]',M.tpart(timeinit));
legendend=sprintf('t=%1.3g [s]',M.tpart(timesteppart));

Tpar=std(0.5*M.me*Vpar.^2*M.vlight^2)/M.qe
Tparend=std(0.5*M.me*Vparend.^2*M.vlight^2)/M.qe
Tperp=std(0.5*M.me*Vperp.^2*M.vlight^2)/M.qe
Tperpend=std(0.5*M.me*Vperpend.^2*M.vlight^2)/M.qe

f=figure;
subplot(1,2,1)
p(1)=scatter(Vpar,Vperp,'.','displayname',legendinit);
hold on
p(2)=scatter(Vparend,Vperpend,'.','displayname',legendend);
xlabel('v_{par}/c')
ylabel('v_\perp/c')
legend('location','southoutside')

zleftlim=1;

vpar=linspace(1,M.vlight,1000);
dN=1e13*M.weight;
dN=0;
rp=9.85e-3;
rm=7e-3;
V=2*pi*(rp^2-rm^2);
N=M.N(:,:,timestepNend);
levels=M.rAthet(Rindex,Zindex)*[1 1];
rAthetposend=find(levels(1)>=M.rAthet(:,zleftlim),1,'last');
rleft=sqrt(levels(1)/(M.rAthet(rAthetposend,1)/M.rgrid(rAthetposend)^2));
rposleft=find(M.rgrid>=rleft,1,'first');
n2=mean(mean(N(rposleft,[1 end],:),3));
vpar2=(dN/V/n2)^2;
if n2==0
    vpar2=0;
end
b=M.rgrid(end);
a=M.rgrid(1);
vd1=((M.potout-M.potinn)*M.phinorm/M.rgrid(Rindex))/M.Bz(Zindex,Rindex)/log(b/a);
vd1=-M.Er(Rindex,Zindex,1)/M.Bz(Zindex,Rindex);
%vd1=-M.Er(Rindex,Zindex,timestepNend)/M.Bz(Zindex,Rindex);
vinit=sqrt(M.kb*10000/M.me);
%vinit=sqrt(M.qe*80/M.me)
%vinit=sqrt(120*M.qe/M.me)
%vd2=((M.potout-M.potinn)*M.phinorm/M.rgrid(25))/M.Bz(1,25)/log(b/a);
[Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,floor(M.nz/2)+1));
Rcurv=griddata(Zmesh,M.rAthet,M.B',M.zgrid(zleftlim),M.rAthet(Rindex,Zindex),'natural')/M.B(Zindex,Rindex);
deltaphicomp=-0.5*M.me/M.qe*(vpar2-(vd1+vinit)^2*(1-M.Rcurv))
ratioparper=vpar2/((vd1+vinit)^2*(1-M.Rcurv))

Ndistrib=M.N(:,:,timestepNend);
model=M.potentialwellmodel(timestepNend);
z=model.z;
r=model.r;
pot=model.pot;
rathet=model.rathet;
Zeval=[M.zgrid(zleftlim) M.zgrid(Zindex)];
Psieval=M.rAthet(Rindex,Zindex)*[1 1];
phis=griddata(Zmesh,M.rAthet,M.pot(:,:,end),Zeval,Psieval,'natural');
deltaphi=-diff(phis)

%deltaphi=-215;
R=griddata(Zmesh,M.rAthet,M.B',M.zgrid(zleftlim),M.rAthet(Rindex,Zindex),'natural')/M.B(Zindex,Rindex);
% if vper is above line then electron is kept axially
vper=sqrt(2*M.qe/M.me*deltaphi/(R-1)+vpar.^2/(R-1))/M.vlight;
vpar=vpar/M.vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
xlimits=xlim;
ylimits=ylim;
if(length(vper)>0)
p(3)=plot(vpar,vper,'b-','displayname','Loss parabolla simul');
hold on
plot(-vpar,vper,'b-')
end

vpar=linspace(1,M.vlight,1000);
vper=sqrt(-2*M.qe/M.me*deltaphicomp/(R-1)+vpar.^2/(R-1))/M.vlight;
vpar=vpar/M.vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
if(length(vper)>0)
p(4)=plot(vpar,vper,'k--','displayname','Loss parabolla prediction');
hold on
plot(-vpar,vper,'k--')
end
xlim(xlimits)
ylim(ylimits)
legend(p)
axis equal

% subplot(2,2,2)
% scatter(Zp,Vpar,'.','displayname','Init')
% hold on
% scatter(Zend,Vparend,'.','displayname','End')
% xlabel('z [m]')
% ylabel('\beta_{par}')
% 
% subplot(2,2,3)
% scatter(Vperp,Rp,'.','displayname','Init')
% hold on
% scatter(Vperpend,Rend,'.','displayname','End')
% xlabel('\beta_\perp')
% ylabel('R [m]')

ax1=subplot(1,2,2);
h=contourf(ax1,M.zgrid,M.rgrid,Ndistrib);
hold on
[r,z]=find(Ndistrib~=0);
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= 'n [m^{-3}]';
Zx=(M.zgrid(Zindex));
Rx=(M.rgrid(Rindex));

plot(ax1,Zx,Rx,'rx','Markersize',12);

sgtitle(sprintf('r=%1.3g [m] z=%1.3g [m] \\Delta\\phi=%1.3g[kV] R=%1.3g',M.rgrid(Rindex),M.zgrid(Zindex),(M.potout-M.potinn)*M.phinorm,M.Rcurv))

M.savegraph(f,sprintf('%s/%s_phasespaceR%dZ%dpos',M.folder,M.name,Rindex,Zindex),[16,12])


%%
f=figure;
p(1)=scatter(Vpar,Vperp,'.','displayname','Initial');
hold on
p(2)=scatter(Vparend,Vperpend,'.','displayname','Final');
xlabel('v_{par}/c')
ylabel('v_\perp/c')


zleftlim=1;

vpar=linspace(1,M.vlight,1000);
dN=1e13*M.weight;
dN=0;
rp=9.85e-3;
rm=7e-3;
V=2*pi*(rp^2-rm^2);
N=M.N(:,:,timestepNend);
levels=M.rAthet(Rindex,Zindex)*[1 1];
rAthetposend=find(levels(1)>=M.rAthet(:,zleftlim),1,'last');
rleft=sqrt(levels(1)/(M.rAthet(rAthetposend,1)/M.rgrid(rAthetposend)^2));
rposleft=find(M.rgrid>=rleft,1,'first');
n2=mean(mean(N(rposleft,[1 end],:),3));
vpar2=(dN/V/n2)^2;
if n2==0
    vpar2=0;
end
b=M.rgrid(end);
a=M.rgrid(1);
vd1=((M.potout-M.potinn)*M.phinorm/M.rgrid(Rindex))/M.Bz(Zindex,Rindex)/log(b/a);
vd1=-M.Er(Rindex,Zindex,1)/M.Bz(Zindex,Rindex);
%vd1=-M.Er(Rindex,Zindex,timestepNend)/M.Bz(Zindex,Rindex);
vinit=sqrt(M.kb*10000/M.me);
%vinit=sqrt(M.qe*80/M.me)
%vinit=sqrt(120*M.qe/M.me)
%vd2=((M.potout-M.potinn)*M.phinorm/M.rgrid(25))/M.Bz(1,25)/log(b/a);
[Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,floor(M.nz/2)+1));

% Find the R magnetic ratio between local position and left
Rcurv=griddata(Zmesh,M.rAthet,M.B',M.zgrid(zleftlim),M.rAthet(Rindex,Zindex),'natural')/M.B(Zindex,Rindex);

deltaphicomp=-0.5*M.me/M.qe*(vpar2-(vd1+vinit)^2*(1-M.Rcurv))

ratioparper=vpar2/((vd1+vinit)^2*(1-M.Rcurv))

Ndistrib=M.N(:,:,timestepNend);
model=M.potentialwellmodel(timestepNend);
z=model.z;
r=model.r;
pot=model.pot;
rathet=model.rathet;
Zeval=[M.zgrid(zleftlim) M.zgrid(Zindex)];
Psieval=M.rAthet(Rindex,Zindex)*[1 1];
phis=griddata(Zmesh,M.rAthet,M.pot(:,:,end),Zeval,Psieval,'natural');
% calculate the potential difference between left and local
deltaphi=-diff(phis)

%deltaphi=-215;
R=griddata(Zmesh,M.rAthet,M.B',M.zgrid(zleftlim),M.rAthet(Rindex,Zindex),'natural')/M.B(Zindex,Rindex);
vper=sqrt(2*M.qe/M.me*deltaphi/(R-1)+vpar.^2/(R-1))/M.vlight;
vpar=vpar/M.vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
xlimits=xlim;
ylimits=ylim;
if(length(vper)>0)
p(3)=plot(vpar,vper,'b-','displayname','Simulation');
hold on
plot(-vpar,vper,'b-')
end

vpar=linspace(1,M.vlight,1000);
vper=sqrt(-2*M.qe/M.me*deltaphicomp/(R-1)+vpar.^2/(R-1))/M.vlight;
vpar=vpar/M.vlight;
vpar=vpar(real(vper)~=0);
vper=vper(real(vper)~=0);
if(length(vper)>0)
p(4)=plot(vpar,vper,'k--','displayname','Prediction');
hold on
plot(-vpar,vper,'k--')
end
xlim(xlimits)
ylim(ylimits)
legend(p)
axis equal
legend('location','northeast')

title(sprintf('r=%1.3g [m] z=%1.3g [m] \\Delta\\phi=%1.3g[kV] R=%1.3g',M.rgrid(Rindex),M.zgrid(Zindex),(M.potout-M.potinn)*M.phinorm,M.Rcurv))

M.savegraph(f,sprintf('%s/%s_phasespaceR%dZ%d',M.folder,M.name,Rindex,Zindex),[12,8])