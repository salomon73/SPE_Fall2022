function PlotParticleTrajectory(M,partnumber,t,text)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if (nargin<4)
    text='';
else
    text=['_',text];
end

if isa(M,'h5parts')
           geomstruct=M.parent;
else
    geomstruct=M;
end

f=figure();
%time=((t(1):t(end))-1)*M.dt*double(M.it2);
time=M.tpart(t(1):t(end));
%sgtitle(sprintf('Particles %d',partnumber));
ax1=subplot(1,2,1);
hold(ax1);
Blines=geomstruct.rAthet;
levels=linspace(min(min(Blines(:,2:end))),max(max(Blines(:,:))),10);
[~,h1]=contour(ax1,geomstruct.zgrid*1000,geomstruct.rgrid*1000,Blines,real(levels),'r-.','linewidth',1.5,'Displayname','Magnetic field lines');

ax2=subplot(1,2,2);
hold(ax2);
Msize=8;
Z=M.Z(partnumber,t,true)*1000;
R=M.R(partnumber,t,true)*1000;
THET=M.THET(partnumber,t,true);
for i=1:length(partnumber)
    Zl=Z(i,R(i,:)>0);
    Rl=R(i,R(i,:)>0);
    THETl=THET(i,R(i,:)>0);
p1(i)=plot(ax1,Zl,Rl,'x-.','Linewidth',1.1,...
    'Displayname',sprintf('part=%d',partnumber(i)),'MarkerIndices',1,...
    'Markersize',Msize);
plot(ax1,Zl(end),Rl(end),'^','Linewidth',1.1,...
    'Displayname',sprintf('part=%d',partnumber(i)),'MarkerIndices',1,...
    'Markersize',Msize,'color',p1(i).Color)
x=Rl.*cos(THETl);
y=Rl.*sin(THETl);
%pp=spline(x,y);
p2=plot(ax2,x,y,'x-.','Linewidth',1.1,...
    'Displayname',sprintf('part=%d',partnumber(i)),'MarkerIndices',1,...
    'Markersize',Msize);
plot(ax2,Rl(end).*cos(THETl(end)),Rl(end).*sin(THETl(end)),'^',...
    'Linewidth',1.1,'Displayname',sprintf('part=%d',partnumber(i)),'MarkerIndices',1,...
    'Markersize',Msize,'color',p2.Color)
end
xlabel(ax1,'z [mm]')
ylabel(ax1,'r [mm]')
sgtitle(sprintf('Position between t=[%3.2f-%3.2f] ns',time(1)*1e9,time(end)*1e9))
xlim(ax1,[M.zgrid(1) M.zgrid(end)]*1e3)
%xlim(ax1,[-.17 .17])
ylim(ax1,[M.rgrid(1) M.rgrid(end)]*1e3)
%ylim(ax1,[0.005 0.025])
grid(ax1,'on');
xlabel(ax2,'x [mm]')
ylabel(ax2,'y [mm]')
t=2*pi*linspace(0,1,100);
if(geomstruct.conformgeom)

xin=geomstruct.rgrid(1)*cos(t)*1e3;
yin=geomstruct.rgrid(1)*sin(t)*1e3;
plot(ax2,xin,yin,'k-')
xout=geomstruct.rgrid(end)*cos(t)*1e3;
yout=geomstruct.rgrid(end)*sin(t)*1e3;
plot(ax2,xout,yout,'k-')
plot(ax1,geomstruct.zgrid*1e3,geomstruct.rgrid(1)*ones(size(geomstruct.zgrid))*1e3,'k-')
plot(ax1,geomstruct.zgrid*1e3,geomstruct.rgrid(end)*ones(size(geomstruct.zgrid))*1e3,'k-')
else
    subplot(1,2,1)
    contour(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3,geomstruct.geomweight(:,:,1),[0 0],'k-')
    [~,rgrid]=meshgrid(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3);
    rlims=rgrid(geomstruct.geomweight(:,:,1)<0);
    xin=geomstruct.rgrid(1)*cos(t)*1e3;
yin=geomstruct.rgrid(1)*sin(t)*1e3;
plot(ax2,xin,yin,'k-')
xout=min(rlims)*cos(t);
yout=min(rlims)*sin(t);
plot(ax2,xout,yout,'k-')
    
end
legend(ax1,p1(1:end))
legend('location','southwest')
axis(ax2,'equal')
%legend(ax2)
grid(ax2,'on');

% ax=subplot(1,2,2);
% xlabel(ax,'t [s]')
% ylabel(ax,'VR [m/s]')
% title(ax,'Radial velocity')
% grid on;

f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.fullpath);
f.PaperUnits='centimeters';
f.PaperSize=[16 9];
print(f,sprintf('%sTrajectories%s',name,text),'-dpdf','-fillpage')
savefig(f,sprintf('%sTrajectories%s',name,text))
end

