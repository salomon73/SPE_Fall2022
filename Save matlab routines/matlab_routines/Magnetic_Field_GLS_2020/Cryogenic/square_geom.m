% creates the magnetic field h5 file necessary for espic2d
% This uses the geometry of the refurbished 170GHz coaxial gyrotron gun 
% The magnetic field is the one created by the asg magnet 
% Both the magnet and geometry are definde in the Final report of the 
% Development of the european gyrotron CCCGDS6 


%% Define the individual boundaries
rmin=0.039;
rmax=0.082;
zmin=0.107;
zmax=0.189;

rc=0.5*(rmin+rmax);
zc=0.5*(zmin+zmax);
zmin=zc-0.5*(rmax-rmin);
zmax=zc+0.5*(rmax-rmin);

npts=9;
geomcells={};
geomcells{1}.Z=zeros(1,4*npts-3);
geomcells{1}.R=zeros(1,4*npts-3);

geomcells{1}.Z(1)=zmax;
geomcells{1}.R(1)=rmin;

i=1;
points=linspace(0,1,npts)'*([zmin,rmin]-[zmax rmin])+[zmax rmin];
geomcells{1}.Z(i+(1:npts-1))=points(2:end,1);
geomcells{1}.R(i+(1:npts-1))=points(2:end,2);

i=i+npts-1;
points=linspace(0,1,npts)'*([zmin,rmax]-[zmin rmin])+[zmin rmin];
geomcells{1}.Z(i+(1:npts-1))=points(2:end,1);
geomcells{1}.R(i+(1:npts-1))=points(2:end,2);

i=i+npts-1;
points=linspace(0,1,npts)'*([zmax,rmax]-[zmin rmax])+[zmin rmax];
geomcells{1}.Z(i+(1:npts-1))=points(2:end,1);
geomcells{1}.R(i+(1:npts-1))=points(2:end,2);

i=i+npts-1;
points=linspace(0,1,npts)'*([zmax,rmin]-[zmax rmax])+[zmax rmax];
geomcells{1}.Z(i+(1:npts-1))=points(2:end,1);
geomcells{1}.R(i+(1:npts-1))=points(2:end,2);





    
geomcells{1}.Z=flip(geomcells{1}.Z);
geomcells{1}.R=flip(geomcells{1}.R);


geomcells{1}.Z=[geomcells{1}.Z(5:end) geomcells{1}.Z(2:5)];
geomcells{1}.R=[geomcells{1}.R(5:end) geomcells{1}.R(2:5)];
for k=1:length(geomcells)
geomcells{k}.order=3;
geomcells{k}.dim=2;
geomcells{k}.epsce=1e-9;
geomcells{k}.epsge=1e-9;
geomcells{k}.Dval=0;
geomcells{k}.name=sprintf('Electrode_%i',k);
geomcells{k}.points=[geomcells{k}.Z; geomcells{k}.R];
geomcells{k}.type=0;
geomcells{k}.periodic=0;
end

%% Plots
f=figure;
for k=1:length(geomcells)
    plothandle=plot(geomcells{k}.Z, geomcells{k}.R,'k-x','linewidth',1.5);
    hold on
    geomcells{k}.points=[geomcells{k}.Z; geomcells{k}.R];
    order=geomcells{k}.order;
    knots=linspace(0,1,length(geomcells{k}.Z)-(order-2));
    knots=augknt(knots, order);
    sizec=size(geomcells{k}.Z);
    order=length(knots)-sizec(end);
    coeffs=[geomcells{k}.Z; geomcells{k}.R];
    pp=spmak(knots,coeffs);
    s=linspace(0,1,1000);
    fittedpos=fnval(pp,s);
    plot(fittedpos(1,:),fittedpos(2,:),'x-')
end
%axis equal
%[~,cont2]=contour(Zphi,Rphi,Phi,20,'b');
%rectangle('Position',[-0.011, 0.06375, 0.032+0.011, 0.081-0.06375],'EdgeColor','magenta','Linestyle','--')

legend([plothandle],{'Gun geometry', 'Magnetic field lines'},'location','southwest')
f.PaperUnits='centimeters';
f.PaperSize=[12,8];
xlabel('z [m]')
ylabel('r [m]')


% print(f,name,'-dpdf','-fillpage')
% savefig(f,name)
% set(f, 'Color', 'w');
% export_fig(f,name,'-eps')
hold off


%% Save magnetic field and geometry to disk
save=true;
overwrite=true;
if save
    savegeomtoh5('test_square_inb.h5',geomcells,1e-2,overwrite);
end

