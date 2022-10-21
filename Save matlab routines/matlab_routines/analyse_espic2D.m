% file='stable_13_decoupled_fft.h5';
file='teststablegauss_13_fft.h5';
% file='teststable_13_1058400.h5';
%close all hidden;
if (~exist('M','var'))
    M.file=file;
end
M=load_espic2d(file,M,'fields');

% M.epot = 0.5*h5read(file,'/data/var0d/epot');
% M.ekin = h5read(file,'/data/var0d/ekin');
% M.etot = M.epot+M.ekin;
% M.eerr = M.etot-M.etot(2);
[~, name, ~] = fileparts(M.file);
tmin=2;
tmax=length(M.ekin);
f=figure();
plot(M.t0d(tmin:tmax),M.ekin(tmin:tmax),'o-',...
    M.t0d(tmin:tmax),M.epot(tmin:tmax),'d-',...
    M.t0d(tmin:tmax),M.etot(tmin:tmax),'h-',...
    M.t0d(tmin:tmax),M.eerr(tmin:tmax),'x--')
legend('E_{kin}', 'E_{pot}', 'E_{tot}','E_{err}')
xlabel('t [s]')
ylabel('Energies [J]')
grid on

f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
f.PaperSize=[15 7];
print(f,sprintf('%s_Energy',name),'-dpdf','-fillpage')

f=figure();
semilogy(M.t0d(tmin:tmax),abs(M.eerr(tmin:tmax)./M.etot(tmin:tmax)),'-')
xlabel('t [s]')
ylabel('|E_{err}/E_{tot}|')
grid on
f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
f.PaperSize=[15 7];
print(f,sprintf('%s_Energyerr',name),'-dpdf','-fillpage')

figure
[Z,R] = meshgrid (M.zgrid,M.rgrid);
subplot(2,1,1)
contourf(Z,R,M.Br')
xlabel('Z')
ylabel('R')
c = colorbar;
c.Label.String= 'B_r [T]';
grid on

subplot(2,1,2)
contourf(Z,R,M.Bz')
xlabel('Z [m]')
ylabel('R [m]')
c = colorbar;
c.Label.String= 'B_z [T]';
grid on

[~, name, ~] = fileparts(M.file);
f=figure('NumberTitle', 'off', 'Name', [name,' Bfield']);
contourf(Z,R,M.B',100,'LineColor','none');
shading flat
colormap('jet')
hold on
Bz=M.Bz';
Br=M.Br';
gridstep=4;
quiver(Z(1:gridstep:end,1:gridstep:end),R(1:gridstep:end,1:gridstep:end),...
       Bz(1:gridstep:end,1:gridstep:end),Br(1:gridstep:end,1:gridstep:end),'k')
% h=streamslice(Z,R,M.Bz',M.Br');
% set(h,'Color','w')
xlabel('Z [m]')
ylabel('R [m]')
c = colorbar;
c.Label.String= '|B| [T]';
f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
f.PaperPosition=[0 0 10 6];
f.PaperSize=[16 10];
print(f,sprintf('%s_Bfield',name),'-dpdf','-fillpage')


[~, name, ~] = fileparts(M.file);
f=figure('NumberTitle', 'off', 'Name', [name,' rAfield']);
subplot(2,1,1)
contour(Z,R,M.Athet'.*R,50);
hold on
xlabel('Z [m]')
ylabel('R [m]')
rectangle('Position',[-.16 4e-3 .32 10e-3],'EdgeColor','r')
subplot(2,1,2)
contour(Z,R,M.Athet'.*R,600);
hold on
xlabel('Z [m]')
ylabel('R [m]')
c = colorbar('SouthOutside');
c.Label.String= '|rA_\theta| [Tm^2]';
ylim([4 14]*10^(-3))
xlim([-0.16 0.16])
f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
f.PaperSize=[14 12];
print(f,sprintf('%s_fieldlines',name),'-dpdf','-fillpage')



% figure
% ax3=subplot(2,1,1);
% contourf(ax3,M.zgrid,M.rgrid,mean(M.Er(:,:,end-2000:end),3)')
% xlabel(ax3,'Z [m]')
% ylabel(ax3,'R [m]')
% c = colorbar(ax3);
% c.Label.String= 'E_r [V/m]';
% title(ax3,'Average Radial Electric field')
% grid on
% ax3.Children(1).ZDataSource='data.Er';
% 
% ax4=subplot(2,1,2);
% contourf(ax4,M.zgrid,M.rgrid,mean(M.Ez(:,:,end-2000:end),3)')
% xlabel(ax4,'Z [m]')
% ylabel(ax4,'R [m]')
% c = colorbar(ax4);
% c.Label.String= 'E_z [V/m]';
% title(ax4,'Average Axial Electric field')
% grid on
% ax4.Children(1).ZDataSource='data.Ez';


f=figure();
tstudied=0;
legtext=sprintf("t=%2.1f - %2.1f [ns]",M.tpart(end-tstudied)*1e9,M.tpart(end)*1e9);
subplot(1,3,1)
VR=M.VR(:,1);
h1=histogram(VR,20,'BinLimits',[min(VR(:)) max(VR(:))],'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
hold on
VR=mean(M.VR(:,end-tstudied:end),2);
h1=histogram(VR,'BinWidth',h1.BinWidth,'DisplayName',legtext);
ylabel('counts')
xlabel('V_r [m/s]')
legend

subplot(1,3,2)
VTHET=M.VTHET(:,1);
h1=histogram(VTHET,20,'BinLimits',[min(VTHET(:)) max(VTHET(:))],'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
hold on
VTHET=mean(M.VTHET(:,end-tstudied:end),2);
h1=histogram(VTHET,'BinWidth',h1.BinWidth,'DisplayName',legtext);
ylabel('counts')
xlabel('V_\theta [m/s]')
legend

subplot(1,3,3)
hold off
VZ=M.VZ(:,1);
h1=histogram(VZ,20,'BinLimits',[min(VZ(:)) max(VZ(:))],'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
hold on
VZ=mean(M.VZ(:,end-tstudied:end),2);
h1=histogram(VZ,'BinWidth',h1.BinWidth,'DisplayName',legtext);
ylabel('counts')
xlabel('V_z [m/s]')
legend
f=gcf;
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);



dispespicFields(M)

BrillouinRatio=2*M.omepe^2/(M.qe*min(min(M.Bz))/M.me)^2