%Mtest=espic2dhdf5('test_fst_80_2d_websurfw.h5',false);
%M=espic2dhdf5('/misc/spltestdwn.h5',false);
Mtest=M;

if length(Mtest.r_a)<1
    r_a=3.e-3;
    r_b=0.152;
    z_0=0;
    r_0=0.158;
    z_r=0.26;
    r_r=0.02;
    Lr=0.12;
    Lz=0.96;
    above=1;
    Interior=-1;
else
    r_a=Mtest.r_a;
    r_b=Mtest.r_b;
    z_0=Mtest.z_0;
    r_0=Mtest.r_0;
    z_r=Mtest.z_r;
    r_r=Mtest.r_r;
    Lr=Mtest.L_r;
    Lz=Mtest.L_z;
    above=double(Mtest.above1);
    Interior=double(Mtest.interior);
end
if isempty(Lr)
Lr=0.12;
Lz=0.96;
end
invr_r=1/r_r;
invr_z=1/z_r;

nbcont=24;

[z,r]=meshgrid(Mtest.zgrid,Mtest.rgrid);

%wgeom=geom_weight(z,r,r_a,above,z_0,r_0,invr_z,invr_r,Interior);
wgeom=Mtest.geomweight(:,:,1);%wgeom(:,:,1);

z0=0.5*(M.zgrid(end)+M.zgrid(1));
r0=0.5*(M.rgrid(end)+M.rgrid(1));

potth=(sin(pi*(z-z0)/Lz).*sin(pi*(r-r0)/Lr)+2);
Erth=-(sin(pi*(z-z0)/Lz).*cos(pi*(r-r0)/Lr)*(pi/Lr));
Ezth=-(cos(pi*(z-z0)/Lz).*sin(pi*(r-r0)/Lr)*(pi/Lz));

% potnum=Mtest.pot(:,:,1)/Mtest.phinorm;
% Er=Mtest.Er(:,:,1)/Mtest.phinorm;
% Ez=Mtest.Ez(:,:,1)/Mtest.phinorm;
potnum=Mtest.potxt/Mtest.phinorm;
Er=Mtest.Erxt/Mtest.phinorm;
Ez=Mtest.Ezxt/Mtest.phinorm;


potnum(wgeom<0)=NaN;
potth(wgeom<0)=NaN;
Erth(wgeom<0)=NaN;
Ezth(wgeom<0)=NaN;
Er(wgeom<0)=NaN;
Ez(wgeom<0)=NaN;
invpotth=1./potth;
%invpotth(isinf(invpotth))=0;
%poterr=abs(poterr./max(potth(:)));
inverth=1./Erth;
inverth(isinf(inverth))=0;
%Ererr=abs(Er-Erth)./max(Erth(:));
%Ezerr=abs(Ez-Ezth)./max(Ezth(:));
Ererr=abs((Er-Erth)./max(Erth(:)));
Ezerr=abs((Ez-Ezth)./max(Ezth(:)));
poterr=abs((potnum-potth)./max(potth(:)));
zerr=z;
rerr=r;
Mtest.name;
L2errpot=norm(potth-potnum,2)/norm(potth,2)
L2errEr=norm(Erth-Er,2)/norm(Erth,2)
L2errEz=norm(Ezth-Ez,2)/norm(Ezth,2)
maxerr=max(max(abs((potth-potnum)./potth)))
gridsize=sprintf(" (nz=%i,nr=%i)",Mtest.nz,sum(Mtest.nnr)); 
  

%% numerical potential solution
f=figure();
lvls=linspace(-1,1,nbcont)*max(abs(potth(:)));
subplot(2,1,1)
contourf(Mtest.zgrid, Mtest.rgrid,potnum,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
caxis([min(potth(:)) max(potth(:))])
title(['pot num',gridsize])
xlabel('z [m]')
ylabel('r [m]')
c = colorbar;
c.Label.String='\phi_{num} [V]';
caxis([min(potth(:)) max(potth(:))])

subplot(2,1,2)
contourf(Mtest.zgrid,Mtest.rgrid,potth,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
c = colorbar;
c.Label.String='\phi_{th} [V]';
caxis([min(potth(:)) max(potth(:))])
title('pot th')
xlabel('z')
ylabel('r')
linkaxes

colormap jet
Mtest.savegraph(f,sprintf('%s/%s_potnum',Mtest.folder,Mtest.name));


%% numerical radial electric field
figure
lvls=linspace(-1,1,nbcont)*max(abs(Erth(:)));
subplot(2,1,1)
contourf(Mtest.zgrid,Mtest.rgrid,Er,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
colorbar
title('E_r')
xlabel('z')
ylabel('r')
subplot(2,1,2)
contourf(Mtest.zgrid,Mtest.rgrid,Erth,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
colorbar
linkaxes
title('E_r th')
xlabel('z')
ylabel('r')

colormap jet

%% numerical axial electric field
figure
lvls=linspace(-1,1,nbcont)*max(abs(Ezth(:)));
subplot(2,1,1)
contourf(Mtest.zgrid,Mtest.rgrid,Ez,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
colorbar
title('E_z')
xlabel('z')
ylabel('r')
subplot(2,1,2)
contourf(Mtest.zgrid,Mtest.rgrid,Ezth,lvls)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
colorbar
linkaxes
title('E_z th')
xlabel('z')
ylabel('r')

colormap jet


%% relative error on potential
f=figure();
poterr(Mtest.geomweight(:,:,1)<0)=0;

contourf(zerr,rerr,poterr,logspace(-8,1,nbcont))
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
set(gca,'ColorScale','log')
%set(gca,'ZScale','log')
title(['Err \phi web',gridsize])
xlabel('z [m]')
ylabel('r [m]')
[poterrmax, id]=max(poterr(2:end-1,2:end-1),[],'all','linear');
[idr,idz]=ind2sub(size(poterr(2:end-1,2:end-1)),id);
fprintf('poterrmax= %1.5f, at z=%1.4f r=%1.4f\n',poterrmax,Mtest.zgrid(idz+1),Mtest.rgrid(idr+1))
c = colorbar;
c.Label.String='relative error';

colormap jet
Mtest.savegraph(f,sprintf('%s/%s_pot_error',Mtest.folder,Mtest.name));

%% relative error on Er
f=figure();
Ererr(Mtest.geomweight(:,:,1)<0)=0;
contourf(zerr,rerr,Ererr,logspace(-8,1,nbcont))
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
set(gca,'ColorScale','log')
%set(gca,'ZScale','log')
title(['Err Er web',gridsize])
xlabel('z [m]')
ylabel('r [m]')
[Ererrmax, id]=max(Ererr(2:end-1,2:end-1),[],'all','linear');
[idr,idz]=ind2sub(size(poterr(2:end-1,2:end-1)),id);
fprintf('Errmax= %1.5f, at z=%1.4f r=%1.4f\n',Ererrmax,Mtest.zgrid(idz+1),Mtest.rgrid(idr+1))
c = colorbar;
c.Label.String='relative error';

colormap jet
Mtest.savegraph(f,sprintf('%s/%s_Er_error',Mtest.folder,Mtest.name));

%% relative error on Ez
f=figure();
Ezerr(Mtest.geomweight(:,:,1)<0)=0;
contourf(zerr,rerr,Ezerr,logspace(-8,1,nbcont))
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
set(gca,'ColorScale','log')
%set(gca,'ZScale','log')
title(['Err Ez web',gridsize])
xlabel('z [m]')
ylabel('r [m]')
[Ezerrmax, id]=max(Ezerr(2:end-1,2:end-1),[],'all','linear');
[idr,idz]=ind2sub(size(poterr(2:end-1,2:end-1)),id);
fprintf('Ezrmax= %1.5f, at z=%1.4f r=%1.4f\n',Ezerrmax,Mtest.zgrid(idz+1),Mtest.rgrid(idr+1))
c = colorbar;
c.Label.String='relative error';
colormap jet
Mtest.savegraph(f,sprintf('%s/%s_Ez_error',Mtest.folder,Mtest.name));



%% analytical potential solution
f=figure();
contourf(Mtest.zgrid, Mtest.rgrid,potth)
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'r-','linewidth',1.2)
title('pot th')
xlabel('z [m]')
ylabel('r [m]')
c = colorbar;
c.Label.String='\phi_{Ana} [V]';
colormap jet
Mtest.savegraph(f,sprintf('%s/%s_potAna',Mtest.folder,Mtest.name));

% %% class of the cells
% zclass=(M.zgrid(1:end-1)+M.zgrid(2:end))/2;
% rclass=(M.rgrid(1:end-1)+M.rgrid(2:end))/2;
% [zclassgrid, rclassgrid]=meshgrid(zclass,rclass);
% gridclass=M.geomweight(:,:,1);%geom_weight(zclassgrid,rclassgrid,r_a,above,z_0,r_0,invr_z,invr_r,Interior);
% 
% f=figure();
% imagesc(zclass,rclass,gridclass(:,:,1))
% hold on
% contour(M.zgrid,M.rgrid,wgeom,[0 0],'r-','linewidth',1.2)

%% geometric weight and boundaries
f=figure();
contourf(Mtest.zgrid,Mtest.rgrid,wgeom,'color','none')
hold on
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'m--','linewidth',2.2)
contour(Mtest.zgrid,Mtest.rgrid,wgeom,max(wgeom(:))*[1 1],'k--','linewidth',2.2)
title(['weight',gridsize])
c = colorbar;
c.Label.String='w(x)';
xlabel('z [m]')
ylabel('r [m]')
%caxis([min(wgeom(:)) 1])
colormap jet
Mtest.savegraph(f,sprintf('%s/%s_weight',Mtest.folder,Mtest.name));


%% cells lines
f=figure();

for i=1:length(Mtest.rgrid)
    plot([Mtest.zgrid(1) Mtest.zgrid(end)],Mtest.rgrid(i)*[1 1],'k-','linewidth',0.3)
    hold on
end
for j=1:length(Mtest.zgrid)
    pl2=plot(Mtest.zgrid(j)*[1 1],[Mtest.rgrid(1) Mtest.rgrid(end)],'k-','linewidth',0.3);
end
contour(Mtest.zgrid,Mtest.rgrid,wgeom,[0 0],'m--','linewidth',2.2)
M.displaysplbound
%legend([f.Children(1).Children(1),pl2],{'Boundaries','Grid lines'},'location','east')
xlabel('z[m]')
ylabel('r[m]')
axis equal

Mtest.savegraph(f,sprintf('%s/%s_geomgrid',Mtest.folder,Mtest.name));

