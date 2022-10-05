% Create a video showing the time evolution of the electron cloud density
% and the equipotential lines
% Draws on top the metallic geometry and the magnetic field lines

% Set the starting and ending time frames
initfieldstep=1;
tend=length(M.t2d);%1000;%length(M.t2d);%min(60000,length(M.t2d));
% set the timestep of drawing
step=2;

% defines the upper limit of the color scale
maxN=2.5e17;
minN=5e12;

%prepare the video writer parameters
filename=[M.name,'_dens_log_presentation.avi'];
videowriterobj=VideoWriter([M.folder,'/',filename]);
videowriterobj.FrameRate=50;
videowriterobj.Quality=90;
open(videowriterobj);

tsize=20;

if M.neutcol.present
    vexb0=(M.Ez(:,:,1).*M.Br'-M.Er(:,:,1).*M.Bz')./(M.B'.^2);
    vexb0(M.geomweight(:,:,1)<=0)=0;
    E=0.5*M.msim/M.weight*mean(abs(vexb0(:)))^2/M.qe;
    taucol=1/(M.neutcol.neutdens*mean(abs(vexb0(:)))*(M.sigio(E)+M.sigmela(E)+M.sigmio(E)));
    try
        Sio_S=1e17*(M.neutcol.neutdens*mean(abs(vexb0(:)))*M.sigio(E))/(M.maxwellsrce.frequency*M.weight/(pi*(M.maxwellsrce.rlim(2)^2-M.maxwellsrce.rlim(1)^2)*diff(M.maxwellsrce.zlim)))
    catch
    end
    tlabel='\tau_d';
else
    taucol=2*pi/M.omece;
    tlabel='\tau_ce';
end




%% Fields
f=figure('Name', sprintf('%s fields',M.name),'Position',[0 0 1600 900]);
ax1=gca;

% extract and plot the initial density
dens=M.N(:,:,initfieldstep);
dispdens=dens;
colormap(flipud(hot));
geomw=M.geomweight(:,:,1);
dispdens(geomw<=0)=0;
Nlvls=logspace( log10(minN), log10(maxN),50);
%Nlvls=linspace(1e12,maxN,50);
[h,curve]=contourf(ax1,M.zgrid*1000,M.rgrid*1000,dispdens,Nlvls,'Displayname','n_e [m^{-3}]');
set(curve,'linecolor','none');
colormap(flipud(hot));
hold on;



% Draws the magnetic field lines
Blines=M.rAthet;
levels=linspace(min(Blines(M.geomweight(:,:,1)>0)),max(Blines(M.geomweight(:,:,1)>0)),20);
Blines(M.geomweight(:,:,1)<0)=NaN;
[~,h1]=contour(ax1,M.zgrid*1000,M.rgrid*1000,Blines,real(levels),'-.','color','k','linewidth',1.5,'Displayname','Magnetic field lines');


% extract and draw the equipotentials
pot=M.pot(:,:,initfieldstep);
pot(geomw<1e-4)=NaN;
pot(geomw<0)=NaN;
potcolor='b';
[c1,h2]=contour(ax1,M.zgrid*1000,M.rgrid*1000,pot/1e3,'--','color',potcolor,'ShowText','on','linewidth',1.2,'Displayname','Equipotentials [kV]');
clabel(c1,h2,'Color',potcolor,'fontsize',tsize,'LabelSpacing',800)

% Draw the metallic boundaries and the geometry itself
[c1,hContour]=contourf(ax1,M.zgrid*1000,M.rgrid*1000,-geomw,[0,0],'linewidth',1.5);
drawnow;
xlim(ax1,[M.zgrid(1)*1000 M.zgrid(end)*1000])
if(M.conformgeom)
    ylim(ax1,[M.rgrid(1)*1000 M.rgrid(rgridend)*1000])
else
    ylim(ax1,[M.rgrid(1)*1000 M.rgrid(end)*1000])
end

% Add the legend and labels
legend([h1,h2],{'Magnetic field lines','Equipotentials [kV]'},'location','southwest')
xlabel(ax1,'z [mm]')
ylabel(ax1,'r [mm]')
c = colorbar(ax1);
c.Label.String= 'Electron density [m^{-3}]';
view(ax1,2)
grid on;
set(ax1,'fontsize',tsize);


% Change the color of the metallic boundaries to grey
hFills=hContour.FacePrims;
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
try
    hFills(end).ColorData = uint8([150;150;150;255]);
    for idx = 1 : numel(hFills)-1
        hFills(idx).ColorData(4) = 0;   % default=255
    end
catch
end

% if necessary extend the display of the inner and outer cylinders
if(M.walltype~=9)
rectangle('Position',[M.zgrid(1)-2e-4 M.rgrid(1)-3e-4 M.zgrid(end)-M.zgrid(1)+4e-4 3e-4]*1000,'FaceColor',[150 150 150]/255,'linewidth',1.5)
rectangle('Position',[M.zgrid(1) M.rgrid(end) M.zgrid(end)-M.zgrid(1) 2e-4]*1000,'FaceColor',[150 150 150]/255,'EdgeColor','none')
ylim([M.rgrid(1)-2e-4 M.rgrid(end)+2e-4]*1000)
end
%rectangle('Position',[M.zgrid(1) M.rgrid(1) 0.06-M.zgrid(1) 0.06-M.rgrid(1)]*1000,'FaceColor',[150 150 150]/255,'EdgeColor','none')

% give the time in the title
title(sprintf('t= %1.2f [ns]',M.t2d(initfieldstep)*1e9))


% if isfield(M.maxwellsrce,'time_start') && M.maxwellsrce.time_start<=M.t2d(1)
%   rectangle('Position',[M.maxwellsrce.zlim(1) M.maxwellsrce.rlim(1) diff(M.maxwellsrce.zlim) diff(M.maxwellsrce.rlim)]*1000,'linestyle','--')
% end

%Adapt the color axis limits and set log scale
set(ax1,'ColorScale','log')
caxis([Nlvls(1) Nlvls(end)]);




% Run the loop to update the electric potential and electron density
for i=initfieldstep:step:tend
    title(ax1,sprintf('t= %1.2f %s',M.t2d(i)/taucol,tlabel),'fontsize',tsize-2)
    dens=M.N(:,:,i);
    dispdens=dens;
    dispdens(geomw<=0)=0;
    pot=M.pot(:,:,i)./1e3;
    pot(geomw<1e-4)=NaN;
    
    
    curve.ZData=dispdens;
    h2.ZData=pot;
    caxis(ax1,[Nlvls(1) Nlvls(end)]);
    %Curve2.YData=pot(rAthetpos,1:end);
    
    
    drawnow;
    hFills=hContour.FacePrims;
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    try
        hFills(1).ColorData = uint8([150;150;150;255]);
        for idx = 2 : numel(hFills)
            hFills(idx).ColorData(4) = 0;   % default=255
        end
    catch
    end
    writeVideo(videowriterobj,getframe(f));
end


close(videowriterobj);