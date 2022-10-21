t2dlength=size(M.t2d,1);
tstep=t2dlength-1;

model=M.potentialwellmodel(tstep);
model2=M.potentialwellmodel(0);
N=M.N(:,:,tstep);
N(M.geomweight(:,:,1)<0)=0;

figure
contourf(M.zgrid,M.rgrid,N);
hold on
contour(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0], 'r')
xlabel('z [m]')
ylabel('r [m]')
[x,y]=ginput(1);
zpos=find(x>M.zgrid,1,'last');
rpos=find(y>M.rgrid,1,'last');
% rpos=68;
% zpos=285;

rAthet_pos=M.rAthet(rpos,zpos);

%%
Psieval=rAthet_pos*ones(M.nz+1,1);
Zeval=M.zgrid;
%Thermal velocity of creation
%vpar0=sqrt(M.kb*22000/M.msim*M.weight);
%fluid
%vpar0=M.fluidUZ(rpos,zpos,t2dlength);

vpar0=sqrt(2/M.me*M.fluidEkin(3,rpos,zpos,tstep));
%% average kinetic energy
%vperp0=sqrt(2/M.me*M.fluidEkin(2,rpos,zpos,tstep));
%% ExB energy
vperp0=(M.Ez(rpos,zpos,tstep)*M.Br(zpos,rpos)-M.Er(rpos,zpos,tstep)*M.Bz(zpos,rpos))/M.B(zpos,rpos)^2;


%rdisp=sort(unique([M.rAthet(:,end);M.rAthet(end,:)]));
[Zmesh,Rmesh]=meshgrid(M.zgrid,M.rgrid);
pot=griddata(model.z,model.rathet,model.pot,Zeval,Psieval);
potxt=griddata(model2.z,model2.rathet,model2.pot,Zeval,Psieval);

[Zinit,~]=meshgrid(M.zgrid,M.rAthet(:,1));
N0=griddata(Zinit,M.rAthet,N,Zeval,Psieval);


% Magnetic field mirror ratio at each grid position
% compared to local position
R=griddata(Zmesh,M.rAthet,M.B',Zeval,Psieval,'natural')/M.B(zpos,rpos);
R=(R-1)*3+1;

% Electrostatic potential on magnetic field line
% coordinates
phis=griddata(Zmesh,M.rAthet,M.pot(:,:,tstep),Zeval,Psieval,'natural')-M.pot(rpos,zpos,tstep);

rposline=griddata(Zmesh,M.rAthet,Rmesh,Zeval,Psieval,'natural');

%%
f=figure('Name','Axial conf');
subplot(2,1,1)
contourf(M.zgrid*1e3,M.rgrid*1e3,N,'edgecolor','none')
hold on
plot(M.zgrid(zpos)*1e3,M.rgrid(rpos)*1e3,'rx')
plot(Zeval*1e3,rposline*1e3,'r--')
contour(M.zgrid*1e3,M.rgrid*1e3,M.geomweight(:,:,1),[0 0],'w-.','linewidth',1.5,'Displayname','Vessel Boundaries');

totpot=M.pot(:,:,tstep);
totpot(M.geomweight(:,:,1)<0)=0;
[c1,h2]=contour(M.zgrid*1e3,M.rgrid*1e3,totpot./1e3,'m--','ShowText','On','linewidth',1.2,'Displayname','Equipotentials [kV]');
%clabel(c1,h2,'Color','white');
clabel(c1,h2,'Color','white');
%ylim([0.077 0.081])
xlabel('z [m]')
ylabel('r [m]')
c=colorbar;
c.Label.String='n [m^{-3}]';
set(gca,'fontsize',12);
drawnow;
texth=h2.TextPrims;
for i=1:length(texth)
    a=texth(i).String;
    y=char(sprintf('%1.2fkV',str2double(a)));
    h2.TextPrims(i).String=y;
end


% subplot(3,1,2)
% plot(M.zgrid*1e3,pot-pot(zpos),'displayname','total well')
% hold on
% plot(M.zgrid*1e3,potxt-potxt(zpos),'displayname','external well')
% ylabel('-e\Delta\phi [eV]')
% yyaxis right
% plot(M.zgrid,N0,'displayname','density')
% ylabel('n [m^{-3}]')
% xlabel('z [m]')
% grid on
% set(gca,'fontsize',12);
% legend('location','northwest')

subplot(2,1,2)
%EparB=0.5*M.msim/M.weight*((M.Ez(rpos,zpos,tstep)*M.Br(zpos,rpos)-M.Er(rpos,zpos,tstep)*M.Bz(zpos,rpos))/M.B(zpos,rpos)^2)^2*(1-R)/M.qe;
EparB=0.5*M.msim/M.weight*(vperp0)^2*(1-R)/M.qe;
Eparphi=-M.qsim/M.weight*phis/M.qe;
Epar=EparB+Eparphi+0.5*M.msim/M.weight/M.qe*vpar0^2;
plot(M.zgrid,Epar,'displayname','E_{par,tot}')
hold on
plot(M.zgrid,EparB,'displayname','\DeltaE_{par,mag}')
plot(M.zgrid,Eparphi,'displayname','\DeltaE_{par,elec}')
plot(M.zgrid,0.5*M.msim/M.weight*(vpar0^2)/M.qe*ones(size(M.zgrid)),':','displayname','E_{par,0}')
grid on
xlabel('z [m]')
ylabel('E_{par} [eV]')
legend('location','north','NumColumns',3)
yyaxis right
plot(M.zgrid,R,'displayname','R')
ylabel('R=B(s)/B(0)')
%ylim(1000*[-1 1])

set(gca,'fontsize',12);

M.savegraph(f,sprintf('%s/%s_Ax_confR%dZ%d_tweaked',M.folder,M.name,rpos,zpos),[12,20])