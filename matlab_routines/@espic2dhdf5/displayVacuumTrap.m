function displayVacuumTrap(obj,time2d)
% Plot the larmor radius for created particles with low energy
% the larmor radius is calculated by considering that the
% initial perpendicular velocity \approx the ExB velocity
% potential well is shown with isolines on top
if iscell(time2d)
    time2d=cell2mat(time2d);
end


[zgrid,rgrid]=ndgrid(obj.zgrid,obj.rgrid);
geomw=obj.geomweight(:,:,1);
inside=griddedInterpolant(zgrid,rgrid,obj.geomweight(:,:,1)');

% Compute the Larmor radius assuming v_perp=v_exb
if time2d>0
    Er=obj.Er(:,:,time2d);
    Ez=obj.Ez(:,:,time2d);
else
    Er=obj.Erxt(:,:,1);
    Ez=obj.Ezxt(:,:,1);
end
rl=abs(obj.me/obj.qe*(-Er.*obj.Bz'+Ez.*obj.Br')./(obj.B.^3)');
% Remove the points outside the geometry
rl(obj.geomweight(:,:,1)<0)=0;

% calculate the normal vector to the magnetic field lines at each position
perpB=zeros(size(obj.sinthet,1),size(obj.sinthet,2),2);
perpB(:,:,1)=sign(-Er).*obj.sinthet; 
perpB(:,:,2)=-sign(-Er).*obj.costhet; 

pos=zeros(size(obj.sinthet,1),size(obj.sinthet,2),2);
pos(:,:,1)=rgrid'; 
pos(:,:,2)=zgrid';

posu=pos+2*rl.*perpB;

isinside=inside(posu(:,:,2),posu(:,:,1))>0;


model=obj.potentialwellmodel(time2d,true);
z=model.z;
r=model.r;
Pot=model.pot;

Potdepth=scatteredInterpolant(z',r',Pot);

effectdepth=0.5*(Potdepth(posu(:,:,2),posu(:,:,1))+Potdepth(pos(:,:,2),pos(:,:,1)));
effectdepth(geomw<0)=NaN;
effectdepth(~isinside)=NaN;
effectdepth(effectdepth<0)=NaN;

f=figure;
contourf(obj.zgrid,obj.rgrid,effectdepth,40,'edgecolor','none')
hold on
contour(obj.zgrid,obj.rgrid,geomw,[0 0],'r-')
N=obj.N(:,:,end);
N=N/max(N(:));
contour(obj.zgrid,obj.rgrid,N,'k--')
colormap(jet)

xlabel('z [m]')
ylabel('r [m]')
c=colorbar;
c.Label.String='Effective well [eV]';

papsize=[14 8];
obj.savegraph(f,sprintf('%s/%s_effwellr_rl_%i',obj.folder,obj.name,floor(mean(time2d))),papsize);


end
