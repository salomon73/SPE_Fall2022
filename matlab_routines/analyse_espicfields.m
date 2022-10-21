%% time extent of plotted variables
% file='teststablegauss_13_traject2.h5';
% file='unigauss.h5';
%  file='Dav_5e14fine_wr+.h5';
% file='stablegauss_8e19fine.h5'
% % %file='teststable_Dav1.h5';
% % %close all hidden;
% if (~exist('M','var'))
%     M.file=file;
%  end
%  M=load_espic2d(file,M,'all');

t2dlength=size(M.t2d,1)
fieldstart=max(1,t2dlength-500);%floor(0.95*t2dlength);
fieldend=t2dlength;%23;
deltat=t2dlength-fieldstart;
%[~, rgridend]=min(abs(M.rgrid-0.045));
rgridend=sum(M.nnr(1:2));
[~, name, ~] = fileparts(M.file);

M.displayenergy



%% Radial profile 
try
    M.displayrprofile([1,fieldstart:fieldend],[],true)
catch
end

%%
%fieldstart=2300; fieldend=2500;
dens=mean(M.N(:,:,fieldstart:fieldend),3);
geomw=M.geomweight(:,:,1);
geomw(geomw<0)=0;
geomw(geomw>0)=max(max(dens));
dispdens=dens;
dispdens(geomw<=0)=NaN;
pot=mean(M.pot(:,:,fieldstart:fieldend),3);
brillratio=2*dens*M.me./(M.eps_0*M.B'.^2);
[R,Z]=meshgrid(M.rgrid,M.zgrid);
Rinv=1./R;
Rinv(isnan(Rinv))=0;
VTHET=mean(M.fluidUTHET(:,:,fieldstart:fieldend),3);
omegare=(VTHET.*Rinv');
maxdens=max(max(dens));


% f=figure();
% ax3=gca;
% surface(ax3,M.zgrid(1:end-1),M.rgrid(1:end-1),(squeeze(mean(M.Presstens(4,:,:,fieldstart:end),4)))*M.vlight)
% xlabel(ax3,'z [m]')
% ylabel(ax3,'r [m]')
% xlim(ax3,[M.zgrid(1) M.zgrid(end)])
% ylim(ax3,[M.rgrid(1) M.rgrid(50)])
% colormap(ax3,'jet')
% c = colorbar(ax3);
% c.Label.String= 'thermal v_\theta [m/s]';
% %c.Limits=[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))];
% %caxis(ax3,[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))])
% title(ax3,'Azimuthal velocity')
% %set(ax3,'colorscale','log')
% grid(ax3, 'on')
% view(ax3,2)
% f.PaperOrientation='landscape';
% [~, name, ~] = fileparts(M.file);
% f.PaperUnits='centimeters';
% papsize=[16 14];
% f.PaperSize=papsize;
% print(f,sprintf('%sfluid_thermtheta',name),'-dpdf','-fillpage')

%% Density
M.displaydensity(fieldstart,fieldend);

%% Fieldswide
f=figure('Name', sprintf('%s fields',name));
ax1=gca;
h=contourf(ax1,M.zgrid,M.rgrid,dispdens,'Displayname','n_e [m^{-3}]');
hold on;
%[c1,h]=contour(ax1,M.zgrid,M.rgrid,pot./1e3,5,'m--','ShowText','off','linewidth',1.5,'Displayname','Equipotentials [kV]');
%clabel(c1,h,'Color','white')
Blines=zeros(size(M.Athet));
for i=1:length(M.zgrid)
    Blines(i,:)=M.Athet(i,:).*M.rgrid';
end
zindex=floor(length(M.zgrid/2));
levels=logspace( log10(min(min(Blines(:,2:end)))), log10(max(max(Blines(:,:)))),8);
contour(ax1,M.zgrid,M.rgrid,Blines',real(levels),'r-.','linewidth',1.5,'Displayname','Magnetic field lines');
[~, hContour]=contourf(ax1,M.zgrid,M.rgrid,geomw,2);
drawnow;
xlim(ax1,[M.zgrid(1) M.zgrid(end)])
ylim(ax1,[M.rgrid(1) M.rgrid(end)])
legend({'n_e [m^{-3}]','Magnetic field lines'},'location','southwest')
xlabel(ax1,'z [m]')
ylabel(ax1,'r [m]')
title(ax1,sprintf('Density t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
c = colorbar(ax1);
c.Label.String= 'n[m^{-3}]';
view(ax1,2)
%set(h,'edgecolor','none');
grid on;
hFills=hContour.FacePrims;
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
hFills(1).ColorData = uint8([255;255;255;255]);
for idx = 2 : numel(hFills)
   hFills(idx).ColorData(4) = 0;   % default=255
end
%rectangle('position',[M.zgrid(5) M.rgrid(3) M.zgrid(end-5)-M.zgrid(5) (M.rgrid(rgridend)-M.rgrid(1))],'Edgecolor','white','linestyle','--','linewidth',2)
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[18 10];
f.PaperSize=papsize;
print(f,sprintf('%sFieldswide',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sFieldswide',name))

%% Fields
f=figure('Name', sprintf('%s fields',name));
ax1=gca;
h=contourf(ax1,M.zgrid*1000,M.rgrid*1000,dispdens,'Displayname','n_e [m^{-3}]');
hold on;
pot(M.geomweight(:,:,1)<0)=0;
Blines=zeros(size(M.Athet));
for i=1:length(M.zgrid)
    Blines(i,:)=M.Athet(i,:).*M.rgrid';
end
zindex=floor(length(M.zgrid/2));
%levels=logspace( log10(min(min(Blines(:,2:end)))), log10(max(max(Blines(:,:)))),10);
levels=linspace(min(min(Blines(:,2:end))),max(max(Blines(:,:))),10);
[~,h1]=contour(ax1,M.zgrid*1000,M.rgrid*1000,Blines',real(levels),'r-.','linewidth',1.5,'Displayname','Magnetic field lines');
max(abs(pot(:)))/1e3
[c1,h2]=contour(ax1,M.zgrid*1000,M.rgrid*1000,pot./1e3,'m--','ShowText','on','linewidth',1.2,'Displayname','Equipotentials [kV]');
clabel(c1,h2,'Color','white')

%contour(ax1,M.zgrid*1000,M.rgrid*1000,M.geomweight(:,:,1),[0 0],'w-','linewidth',1.5);

[c1,hContour]=contourf(ax1,M.zgrid*1000,M.rgrid*1000,geomw);
drawnow;
xlim(ax1,[M.zgrid(1)*1000 M.zgrid(end)*1000])
if(M.conformgeom)
    ylim(ax1,[M.rgrid(1)*1000 M.rgrid(rgridend)*1000])
else
    ylim(ax1,[M.rgrid(1)*1000 M.rgrid(end)*1000])
end
legend([h1,h2],{'Magnetic field lines','Equipotentials [kV]'},'location','southwest')
xlabel(ax1,'z [mm]')
ylabel(ax1,'r [mm]')
%title(ax1,sprintf('Density t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
c = colorbar(ax1);
c.Label.String= 'n[m^{-3}]';
view(ax1,2)
%set(h,'edgecolor','none');
grid on;
hFills=hContour.FacePrims;
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
try
hFills(1).ColorData = uint8([150;150;150;255]);
for idx = 2 : numel(hFills)
   hFills(idx).ColorData(4) = 0;   % default=255
end
catch
end
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
if M.maxwellsrce.present
    rlen=diff(M.maxwellsrce.rlim);
    zlen=diff(M.maxwellsrce.zlim);
rectangle('Position',[M.maxwellsrce.zlim(1) M.maxwellsrce.rlim(1) zlen rlen]*1000,'Edgecolor','g','Linewidth',3,'Linestyle','--')
end
%rectangle('Position',[-10 73 44 8],'Edgecolor','g','Linewidth',3,'Linestyle','--')
f.PaperUnits='centimeters';
papsize=[12 8];
f.PaperSize=papsize;
rectangle('Position',[M.zgrid(1) M.rgrid(1) M.zgrid(end)-M.zgrid(1) 0.064-M.rgrid(1)]*1000,'FaceColor',[150 150 150]/255)
print(f,sprintf('%sFields',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sFields',name))

%% Brillouin
f=figure('Name', sprintf('%s Brillouin',name));
ax1=gca;
brillratio(geomw<=0)=NaN;
h=surface(ax1,M.zgrid,M.rgrid,brillratio);
set(h,'edgecolor','none');
xlim(ax1,[M.zgrid(1) M.zgrid(end)])
if(M.conformgeom)
    ylim(ax1,[M.rgrid(1) M.rgrid(rgridend)])
else
    ylim(ax1,[M.rgrid(1) M.rgrid(end)])
end
xlabel(ax1,'z [m]')
ylabel(ax1,'r [m]')
title(ax1,'Position')
c = colorbar(ax1);
c.Label.String= 'Brillouin Ratio';
view(ax1,2)
caxis([0 max(1.2,max(brillratio(:)))])
title(ax1,sprintf('Brillouin Ratio t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%sfluid_Brillouin',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sfluid_Brillouin',name))

%% Omegar
f=figure('Name', sprintf('%s Omegar',name));
ax3=gca;
omegare(geomw<=0)=NaN;
surface(ax3,M.zgrid,M.rgrid,omegare,'edgecolor','None')
xlabel(ax3,'z [m]')
ylabel(ax3,'r [m]')
xlim(ax3,[M.zgrid(1) M.zgrid(end)])
if(M.conformgeom)
    ylim(ax3,[M.rgrid(1) M.rgrid(rgridend)])
else
    ylim(ax3,[M.rgrid(1) M.rgrid(end)])
end
colormap(ax3,'jet')
c = colorbar(ax3);
c.Label.String= '|\omega_\theta| [1/s]';
%c.Limits=[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))];
%caxis(ax3,[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))])
title(ax3,sprintf('Azimuthal frequency t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
%set(ax3,'colorscale','log')
grid(ax3, 'on')
view(ax3,2)
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%sfluid_omegar',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sfluid_omegar',name))

%% VTHET
f=figure('Name', sprintf('%s Vthet',name));
ax3=gca;
VTHET(geomw<=0)=NaN;
surface(ax3,M.zgrid,M.rgrid,VTHET,'edgecolor','None')
xlabel(ax3,'z [m]')
ylabel(ax3,'r [m]')
xlim(ax3,[M.zgrid(1) M.zgrid(end)])
if(M.conformgeom)
    ylim(ax3,[M.rgrid(1) M.rgrid(rgridend)])
else
    ylim(ax3,[M.rgrid(1) M.rgrid(end)])
end
colormap(ax3,'jet')
c = colorbar(ax3);
c.Label.String= '|V_\theta| [m/s]';
%c.Limits=[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))];
%caxis(ax3,[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))])
title(ax3,sprintf('Azimuthal velocity t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
%set(ax3,'colorscale','log')
grid(ax3, 'on')
view(ax3,2)
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%sfluid_vthet',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sfluid_vthet',name))

%% vdrift
f=figure('Name', sprintf('%s V drift',name));
ax3=gca;
vdrift=(mean(M.Ez(:,:,fieldstart:end),3).*M.Br'-mean(M.Er(:,:,fieldstart:end),3).*M.Bz')./(M.B.^2)';
vdrift(geomw<=0)=NaN;
v=vdrift;
surface(ax3,M.zgrid,M.rgrid,v,'edgecolor','None')
xlabel(ax3,'z [m]')
ylabel(ax3,'r [m]')
xlim(ax3,[M.zgrid(1) M.zgrid(end)])
if(M.conformgeom)
    ylim(ax3,[M.rgrid(1) M.rgrid(rgridend)])
else
    ylim(ax3,[M.rgrid(1) M.rgrid(end)])
end
colormap(ax3,'jet')
c = colorbar(ax3);
set(ax3,'zscale','log')
set(ax3,'colorscale','log')
c.Label.String= 'V_{ExB}';
%c.Limits=[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))];
%caxis(ax3,[min(M.fluidUTHET(:)) max(M.fluidUTHET(:))])
title(ax3,sprintf('V_{ExB} velocity t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
%set(ax3,'colorscale','log')
grid(ax3, 'on')
view(ax3,2)
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%sfluid_EB',name),'-dpdf','-fillpage')
savefig(f,sprintf('%sfluid_EB',name))


%% HP
if(M.R.nt>=2)
taveragestart=max(floor(0.8*size(M.tpart,1)),2);
M.displayHP(taveragestart);
end

%% VR, VTHET, VZ
M.displayVdistribRThetZ()

%% Vperp Vpar
M.displayVdistribParPer()


%% 0D diagnostics
BrillouinRatioth=2*M.omepe^2/M.omece^2
BrillouinMax=max(max(brillratio))
NpartsTrapped=mean(M.npart(end-10:end))*M.msim/M.me
dblengthiter=fieldstart:fieldend;


Pr=squeeze(M.Presstens(1,:,:,dblengthiter));
Pz=squeeze(M.Presstens(6,:,:,dblengthiter));
enddens=M.N(:,:,dblengthiter);
meanPr=mean(Pr,3);
meanPz=mean(Pz,3);
Tr=Pr./enddens;
Tz=Pz./enddens;
Debye_length=sqrt(abs(min(min(min(mean(Tr(Tr>0)./enddens(Tr>0),3))),min(min(mean(Tz(Tz>0)./enddens(Tz>0),3)))))*M.eps_0/M.qe^2)


%% Debye length
f=figure('Name', sprintf('%s Debye',name));
subplot(2,1,1)
surface(M.zgrid,M.rgrid,sqrt(mean(abs(Tr./enddens*M.eps_0/M.qe^2),3)),'edgecolor','none')
colorbar
f.PaperOrientation='landscape';
xlabel('z [m]')
ylabel('r [m]')
c = colorbar;
c.Label.String= '\lambda_r [m]';

subplot(2,1,2)
surface(M.zgrid,M.rgrid,sqrt(mean(abs(Tz./enddens*M.eps_0/M.qe^2),3)),'edgecolor','none')
colorbar
f.PaperOrientation='landscape';
xlabel('z [m]')
ylabel('r [m]')
c = colorbar;
c.Label.String= '\lambda_z [m]';

[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%s_dblength',name),'-dpdf','-fillpage')
savefig(f,sprintf('%s_dblength',name))

%% Temperature
f=figure('Name', sprintf('%s Temperature',name));
ax1=subplot(2,1,1);
title('Radial temperature')
surface(M.zgrid,M.rgrid,mean(Tr,3)/M.qe,'edgecolor','none')
colormap('jet')
c = colorbar;
c.Label.String= 'T_{er} [eV]';
templimits1=caxis;

ax2=subplot(2,1,2);
title('Axial tempreature')
surface(M.zgrid,M.rgrid,mean(Tz,3)/M.qe,'edgecolor','none')
colormap('jet')
c = colorbar;
c.Label.String= 'T_{ez} [eV]';
maxtemp=max([templimits1,caxis]);
caxis(ax1,[0 maxtemp])
caxis(ax2,[0 maxtemp])
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%s_temp',name),'-dpdf','-fillpage')
savefig(f,sprintf('%s_temp',name))


%% Wells
M.display2Dpotentialwell(t2dlength-1);
M.display2Dpotentialwell(t2dlength-1,false);



