

z=linspace(-0.05,0.2,100);
r=linspace(0.02,0.09,100);
[Z,R] = meshgrid(z,r);

% I = [ 9.7 7.7 96.2499 91.0124 87.4723]; 
% I = [0  0  6.79541  1.967  88.2799];
% I = [0  0  5  1.967  88.2799];
I = [61.56  66.45  113.43  108.09 ]; 
res=zeros([size(Z),3]);
%Bz=zeros(size(Z));
%Br=zeros(size(Z));

for i=1:size(R,1)
res(i,:,:) = B_Ellip_Cryogenic_170({'aphi','bz','br'},'cryogenic',I,R(i,:),Z(i,:));
%Bz(i,:)   = B_Ellip_Cryogenic_170('bz','asg_modified_201206',I,R(i,:),Z(i,:));
%Br(i,:)   = B_Ellip_Cryogenic_170('br','asg_modified_201206',I,R(i,:),Z(i,:));
end
Aphi=res(:,:,1);
Bz=res(:,:,2);
Br=res(:,:,3);


% zphi=linspace(0.005,0.0195,100);
% rphi=linspace(0.06375,0.079,150);
% b=max(rphi);
% a=min(rphi);
% phib=0;
% phia=-30000;
% phi=((phib-phia)*log(rphi)+phia*log(b)-phib*log(a))/log(b/a);
% Er=((phib-phia)./rphi)/log(b/a);
% [Zphi,Rphi] = meshgrid(zphi,rphi);
% Phi=repmat(phi',1,length(zphi));


%% Define the individual boundaries
geom=importproto('../geometry_proto.data',[2,Inf]);
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
    geomcells{j}.Z(n)=geom.Z(i);
    geomcells{j}.R(n)=geom.R(i);
    i=i+1;
    n=n+1;
end
geomcells{1}.name='Cathode';
geomcells{1}.Dval=-30000;
geomcells{2}.name='Body';
geomcells{2}.Dval=0;
geomcells{3}.name='Coaxial insert';
geomcells{3}.Dval=0;
geomcells{4}.name='insulator left';
geomcells{4}.Dval=0;
geomcells{5}.name='insulator right';
geomcells{5}.Dval=0;
for k=1:length(geomcells)
geomcells{k}.order=3;
geomcells{k}.dim=2;
geomcells{k}.epsce=1e-9;
geomcells{k}.epsge=1e-9;
geomcells{k}.type=0;
geomcells{k}.periodic=0;
end
geomcells{end-1}.type=2;
geomcells{end}.type=2;

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
    s=linspace(0,1,10000);
    fittedpos=fnval(pp,s);
    plot(fittedpos(1,:),fittedpos(2,:),'x-')
end
hold on
%end
%axis equal
[~,cont1]=contour(Z,R,R.*Aphi,15,'r:','linewidth',3);
%xlim([min(z) 0.04])
%ylim([0.055 0.085])
%[~,cont2]=contour(Zphi,Rphi,Phi,20,'b');
%rectangle('Position',[-0.011, 0.06375, 0.032+0.011, 0.081-0.06375],'EdgeColor','magenta','Linestyle','--')



legend([plothandle,cont1],{'Gun geometry', 'Magnetic field lines','Wall approximation'},'location','southwest')
f.PaperUnits='centimeters';
f.PaperSize=[12,8];
xlabel('z [m]')
ylabel('r [m]')
print(f,'phiBprofile_protoasg','-dpdf','-fillpage')
savefig(f,'phiBprofile_protoasg')
hold off


%% Save magnetic field and geometry to disk
save=true;
overwrite=true;
if save
%     idr=1:1:length(magnet.r);
% idz=1:1:length(magnet.z);
% Aphi=zeros(length(idr),length(idz));
% Bz=zeros(length(idr),length(idz));
% Br=zeros(length(idr),length(idz));
% 
% for i=1:size(magnet.subcoils,1)
%     Aphi=Aphi+Icoil(i)*magnet.subcoils{i,1}(idr,idz);
%     Bz=Bz+Icoil(i)*magnet.subcoils{i,2}(idr,idz);
%     Br=Br+Icoil(i)*magnet.subcoils{i,3}(idr,idz);
% end
%     savemagtoh5([name,'.h5'],magnet.r,magnet.z,Aphi,Br,Bz,overwrite);
    savegeomtoh5('proto_geom.h5',geomcells,1e-2,overwrite);
end
