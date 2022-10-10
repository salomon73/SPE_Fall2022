function displayLarmorRad2Dwell(obj,time2d)
% Plot the larmor radius for created particles with low energy
% the larmor radius is calculated by considering that the
% initial perpendicular velocity \approx the ExB velocity
% potential well is shown with isolines on top
if iscell(time2d)
    time2d=cell2mat(time2d);
end
if nargin<3
    clims=[-inf inf];
end
f=figure('Name',sprintf('%s Potential well + Larmor radius',obj.name));
ax1=gca;
model=obj.potentialwellmodel(time2d);
z=model.z;
r=model.r;
Pot=model.pot;
rathet=model.rathet;
if (time2d==0)
    title(sprintf('\\rho_l and Potential well Vacuum'))
else
    title(sprintf('\\rho_l and Potential well t=%1.2f [ns]',obj.t2d(time2d)*1e9))
end


geomw=obj.geomweight(:,:,1);


if time2d>0
    Er=obj.Er(:,:,time2d);
    Ez=obj.Ez(:,:,time2d);
else
    Er=obj.Erxt(:,:,1);
    Ez=obj.Ezxt(:,:,1);
end

rl=abs(obj.me/obj.qe*(-Er.*obj.Bz'+Ez.*obj.Br')./(obj.B.^3)');
rl(obj.geomweight(:,:,1)<0)=0;
contourf(obj.zgrid*1e3,obj.rgrid*1e3,rl*1e3)
hold on
%contour(obj.zgrid*1e3,obj.rgrid*1e3,obj.geomweight(:,:,1),[0 0],'r-','linewidth',3)
c=colorbar;
c.Label.String='r_L [mm]';

id=find(time2d==0);
time2d(id)=[];

[Zmesh,Rmesh]=meshgrid(obj.zgrid,obj.rgrid);
Pot=griddata(z,r,Pot,Zmesh,Rmesh,'natural');
%Pot(obj.geomweight(:,:,1)<=0)=NaN;
Pot=Pot/max(max(Pot(~isnan(Pot))))*max(rl(:))*1e3;
contour(obj.zgrid(1:end)*1e3,obj.rgrid(1:end)*1e3,Pot(1:end,1:end),10,'r--','Displayname','Well')
xlabel('z [mm]')
ylabel('r [mm]')
xlim([obj.zgrid(1) obj.zgrid(end)]*1e3)
ylim([obj.rgrid(1) obj.rgrid(end)]*1e3)
hold(gca, 'on')

rdisp=obj.rgrid;

%% Magnetic field lines
Blines=obj.rAthet;
levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
Blines(obj.geomweight(:,:,1)<0)=NaN;
[~,h1]=contour(obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'k-.','linewidth',0.8,'Displayname','Magnetic field lines');

c=colorbar;
colormap('jet');

% Grey outline showing the metalic walls
geomw(obj.geomweight(:,:,1)>0)=-1;
geomw(obj.geomweight(:,:,1)<=0)=1;
[c1,hContour]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,geomw, [0 0]);

drawnow;
xlim(ax1,[obj.zgrid(1)*1000 obj.zgrid(end)*1000])
if(obj.conformgeom)
    ylim([ax1 ],[obj.rgrid(1)*1000 obj.rgrid(rgridend)*1000])
else
    ylim([ax1],[obj.rgrid(1)*1000 obj.rgrid(end)*1000])
end
%ylim(ax1,[0.05*1000 obj.rgrid(end)*1000])
%xlim([obj.zgrid(1) 0.185]*1e3)
xlabel(ax1,'z [mm]')
ylabel(ax1,'r [mm]')
view(ax1,2)
c.Label.String='\rho_L [mm]';
f.PaperUnits='centimeters';
caxis(clims)
grid on;
hFills=hContour.FacePrims;
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
try
    drawnow
    hFills(1).ColorData = uint8([150;150;150;255]);
    for idx = 2 : numel(hFills)
        hFills(idx).ColorData(4) = 0;   % default=255
    end
catch
end

% add central and external metallic walls if we have a coaxial
% configuration
if( obj.walltype >=2 && obj.walltype<=4)
    rectangle('Position',[obj.zgrid(1) obj.r_b obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
    ylimits=ylim;
    ylim([ylimits(1),ylimits(2)+1])
end
if sum(obj.geomweight(:,1,1))==0
    rectangle('Position',[obj.zgrid(1) obj.r_a-0.001 obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
    ylimits=ylim;
    ylim([ylimits(1)-1,ylimits(2)])
end

%axis equal
%xlim([-100 200])
[max_depth,id]=max(abs(Pot(:)));
[idr,idz]=ind2sub(size(Pot),id);
fprintf('Maximum potential wel depth: %f eV\n',max_depth)
fprintf('at location r=%f z=%f [mm]\n',obj.rgrid(idr)*1e3, obj.zgrid(idz)*1e3)
papsize=[14 8];
obj.savegraph(f,sprintf('%s/%s_wellr_rl_%i',obj.folder,obj.name,floor(mean(time2d))),papsize);


end
