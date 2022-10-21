% creates the magnetic field h5 file necessary for espic2d
% This uses the geometry of the refurbished 170GHz coaxial gyrotron gun 
% The magnetic field is the one created by the asg magnet 
% Both the magnet and geometry are definde in the Final report of the 
% Development of the european gyrotron CCCGDS6 


%% Define the individual boundaries
geom=importrefurb('../test_slanted_geom.data');
geomcells={};
j=1;
i=1;
n=1;
while i<=size(geom.Z,1)
    if isnan(geom.Z(i))
        j=j+1;
        while i<=size(geom.Z,1) && isnan(geom.Z(i))
            i=i+1;
        end
        n=1;
    end
    if i<=size(geom.Z,1)
    geomcells{j}.Z(n)=geom.Z(i)/1e3;
    geomcells{j}.R(n)=geom.R(i)/1e3;
    i=i+1;
    n=n+1;
    end
end
for k=1:length(geomcells)
geomcells{k}.order=4;
geomcells{k}.dim=2;
geomcells{k}.epsce=1e-9;
geomcells{k}.epsge=1e-9;
geomcells{k}.Dval=0;
geomcells{k}.name=sprintf('Electrode_%i',k);
geomcells{k}.points=[geomcells{k}.Z; geomcells{k}.R];
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
%geomcellsorig{1}.name='Cathode';
%geomcellsorig{1}.Dval=0;
for k=1:length(geomcellsorig)
geomcellsorig{k}.order=3;
geomcellsorig{k}.dim=2;
geomcellsorig{k}.epsce=1e-9;
geomcellsorig{k}.epsge=1e-9;
end

%% Plots
f=figure;
for k=1:length(geomcells)
    plothandle=plot(geomcells{k}.Z, geomcells{k}.R,'k-','linewidth',1.5);
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
for k=1:length(geomcellsorig)
    plothandleorig=plot(geomcellsorig{k}.Z, geomcellsorig{k}.R,'r--','linewidth',1.5);
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
    savegeomtoh5('test_slanted_geom_sep.h5',geomcells,1e-2,overwrite);
end

