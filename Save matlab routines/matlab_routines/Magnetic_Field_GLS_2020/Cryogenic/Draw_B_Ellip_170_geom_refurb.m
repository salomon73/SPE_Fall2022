% creates the magnetic field h5 file necessary for espic2d
% This uses the geometry of the refurbished 170GHz coaxial gyrotron gun 
% The magnetic field is the one created by the asg magnet 
% Both the magnet and geometry are definde in the Final report of the 
% Development of the european gyrotron CCCGDS6 

%% Import and calculate the magnetic field
magnet=load('asg_magnet_red.mat');



%% calculate magnetic field and magnetic vector
idr=1:10:length(magnet.r);
idz=1:10:length(magnet.z);

% Define the coils currents
%I = [ 9.7 7.7 96.2499 91.0124 87.4723];
% I = [0  0  4.5  1.967  88.2799];
% name='phiBprofile_refurbasg_4_5_red';
% I = [0  0  4.8  1.967  88.2799];
% name='phiBprofile_refurbasg_well_red';
% I = [0  0  5  1.967  88.2799];
% name='phiBprofile_refurbasg_5_red';
% I = [0  0  5.15  1.967  88.2799];
% name='phiBprofile_refurbasg_limdwn_red';
%I = [0  0  5.35  1.967  88.2799];
%name='phiBprofile_refurbasg_limup_red';
% I = [0  0  5.5  1.967  88.2799];
% name='phiBprofile_refurbasg_5_5_red';
%I = [0  0  5.6  1.967  88.2799];
%name='phiBprofile_refurbasg_5_6_red';
%I = [0  0  5.7  1.967  88.2799];
%name='phiBprofile_refurbasg_5_7_red';
% I = [0  0  5.75  1.967  88.2799];
% name='phiBprofile_refurbasg_5_75_red';
% I = [0  0  5.8  1.967  88.2799];
% name='phiBprofile_refurbasg_5_8_red';
% I = [0  0  5.85  1.967  88.2799];
% name='phiBprofile_refurbasg_5_85_red';
% I = [0  0  5.9  1.967  88.2799];
% name='phiBprofile_refurbasg_5_9_red';
% I = [0  0  6  1.967  88.2799];
% name='phiBprofile_refurbasg_6_red';

% I = [0  0  6.4  1.967  88.2799];
% name='phiBprofile_refurbasg_6_4_red';
% I = [0  0  6.79541  1.967  88.2799];
% name='phiBprofile_refurbasg_nominal_red';

% for Ib=linspace(4.6,6.8,12)

Ib=4.4;
I = [0  0  Ib  1.967  88.2799];
name=sprintf('phiBprofile_refurbasg_%2d_red',floor(Ib*10));

Icoil=[I(1)  I(2) -I(3)-I(5) -I(3)-I(5) I(4)+I(5) I(4)+I(5) I(5) I(5)];
rmin    = [0.13958    0.13958    0.13710    0.16952    0.13670    0.19101   0.13458   0.17197   ];
rmax    = [0.142694   0.142694   0.16733    0.20005    0.18941    0.23666   0.16972   0.20412   ];
zmin    = [0.07854    0.13854    0.18200    0.18200    0.31778    0.31748   0.50685   0.50685   ];
zmax    = [0.09846    0.15846    0.28350    0.28350    0.48977    0.49017   0.71128   0.71128   ];
% Nturn   = [    218        218     2528.5       2626       7689    11517.5      6110    9656.5   ];
% na      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils (radial direction
% nb      = 2. *[  10         10         10          10          10          10         10         10   ];   % Number of subcoils in longitudinal direction
% Icoil=Icoil.*Nturn./(na.*nb);

[Z,R] = meshgrid(magnet.z(idz),magnet.r(idr));

Aphi=zeros(length(idr),length(idz));
Bz=zeros(length(idr),length(idz));
Br=zeros(length(idr),length(idz));

for i=1:size(magnet.subcoils,1)
    Aphi=Aphi+Icoil(i)*magnet.subcoils{i,1}(idr,idz);
    Bz=Bz+Icoil(i)*magnet.subcoils{i,2}(idr,idz);
    Br=Br+Icoil(i)*magnet.subcoils{i,3}(idr,idz);
end

%% Define the individual boundaries
geom=importrefurb('../refurb_modif.data');
geomcells={};
j=1;
i=1;
n=1;
while i<=size(geom.Z,1)
    if isnan(geom.Z(i))
        j=j+1;
        while isnan(geom.Z(i))
            i=i+1;
        end
        n=1;
    end
    geomcells{j}.Z(n)=geom.Z(i)/1e3;
    geomcells{j}.R(n)=geom.R(i)/1e3;
    i=i+1;
    n=n+1;
end
geomcells{1}.name='Cathode';
geomcells{1}.Dval=-30000;
geomcells{2}.name='Body';
geomcells{2}.Dval=0;
geomcells{3}.name='Coaxial insert';
geomcells{3}.Dval=0;
for k=1:length(geomcells)
geomcells{k}.order=3;
geomcells{k}.dim=2;
geomcells{k}.epsce=1e-9;
geomcells{k}.epsge=1e-9;
geomcells{k}.type=0;
geomcells{k}.periodic=0;
end


%% Load the original boundary
geomorig=importrefurb('../refurb.data');
geomcellsorig={};
j=1;
i=1;
n=1;
while i<=size(geomorig.Z,1)
    if isnan(geomorig.Z(i))
        j=j+1;
        while isnan(geomorig.Z(i))
            i=i+1;
        end
        n=1;
    end
    geomcellsorig{j}.Z(n)=geomorig.Z(i)/1e3;
    geomcellsorig{j}.R(n)=geomorig.R(i)/1e3;
    i=i+1;
    n=n+1;
end
geomcellsorig{1}.name='Cathode';
geomcellsorig{1}.Dval=-30000;
geomcellsorig{2}.name='Body';
geomcellsorig{2}.Dval=0;
geomcellsorig{3}.name='Coaxial insert';
geomcellsorig{3}.Dval=0;
for k=1:length(geomcellsorig)
geomcellsorig{k}.order=3;
geomcellsorig{k}.dim=2;
geomcellsorig{k}.epsce=1e-9;
geomcellsorig{k}.epsge=1e-9;
end

%% Plots
f=figure;
for k=1:length(geomcellsorig)
    plothandleorig=plot(geomcellsorig{k}.Z*1e3, geomcellsorig{k}.R*1e3,'r-','linewidth',2);
    hold on
end
for k=1:length(geomcells)
    %plothandle=plot(geomcells{k}.Z, geomcells{k}.R,'k-','linewidth',1.5);
    hold on
    geomcells{k}.points=[geomcells{k}.Z; geomcells{k}.R];
    order=geomcells{k}.order;
    knots=linspace(0,1,length(geomcells{k}.Z)-(order-2));
    knots=augknt(knots, order);
    sizec=size(geomcells{k}.Z);
    order=length(knots)-sizec(end);
    coeffs=[geomcells{k}.Z; geomcells{k}.R];
    pp=spmak(knots,coeffs);
    s=linspace(0,1,500);
    fittedpos=fnval(pp,s);
    plot(fittedpos(1,:)*1e3,fittedpos(2,:)*1e3,'--','linewidth',2.2)
end

%axis equal
raphi=R.*Aphi;
lvls=logspace(-4,log10(max(raphi(:))),50);
%[~,cont1]=contour(Z*1e3,R*1e3,raphi,lvls,'b:','linewidth',1.5);
%xlim([min(magnet.z) max(magnet.z)]*1e3)
%ylim([min(magnet.r) max(magnet.r)]*1e3)
%[~,cont2]=contour(Zphi,Rphi,Phi,20,'b');
rectangle('Position',[-0.1, 0.045, 0.292, 0.092-0.045]*1e3,'EdgeColor','black','Linestyle','--','linewidth',2)

%legend([plothandle,cont1],{'Gun geometry', 'Magnetic field lines'},'location','southwest')
%f.PaperUnits='centimeters';
%f.PaperSize=[12,8];
xlabel('z [mm]')
ylabel('r [mm]')
axis equal
grid on
ylim([0,inf])
xlim([-300 750])
for i=1:length(rmin)
    rectangle('Position',[zmin(i) rmin(i) zmax(i)-zmin(i) rmax(i)-rmin(i)]*1e3,'edgecolor','r')
    text(zmin(i)*1e3-10,(rmin(i)+rmax(i))/2*1e3,sprintf('S%i',i),'fontsize',10,'color','r')
end

print(f,name,'-dpdf','-fillpage')
savefig(f,name)
set(f, 'Color', 'w');
export_fig(f,name,'-eps')
hold off



%% Save magnetic field and geometry to disk
save=true;
overwrite=true;
if save
    idr=1:1:length(magnet.r);
idz=1:1:length(magnet.z);
Aphi=zeros(length(idr),length(idz));
Bz=zeros(length(idr),length(idz));
Br=zeros(length(idr),length(idz));

for i=1:size(magnet.subcoils,1)
    Aphi=Aphi+Icoil(i)*magnet.subcoils{i,1}(idr,idz);
    Bz=Bz+Icoil(i)*magnet.subcoils{i,2}(idr,idz);
    Br=Br+Icoil(i)*magnet.subcoils{i,3}(idr,idz);
end
    savemagtoh5([name,'.h5'],magnet.r,magnet.z,Aphi,Br,Bz,overwrite);
    %savegeomtoh5('refurb_geom.h5',geomcells,1e-2,overwrite);
end
% end

