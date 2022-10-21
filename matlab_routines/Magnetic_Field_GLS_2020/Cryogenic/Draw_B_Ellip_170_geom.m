geom=importgeometry('../geometry.data');

% z=linspace(-0.012,0.04,100);
% r=linspace(0.063,0.085,50);
% [Z,R] = meshgrid(z,r);
% 
% I = [61.56  66.45  113.43  108.09 ]; 
% Aphi=zeros(size(Z));
% Bz=zeros(size(Z));
% Br=zeros(size(Z));
% 
% for i=1:size(R,1)
% Aphi(i,:) = B_Ellip_Cryogenic_170('aphi','cryogenic',I,R(i,:),Z(i,:));
% Bz(i,:)   = B_Ellip_Cryogenic_170('bz','cryogenic',I,R(i,:),Z(i,:));
% Br(i,:)   = B_Ellip_Cryogenic_170('br','cryogenic',I,R(i,:),Z(i,:));
% end

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



f=figure;
for k=1:size(geom.Z)
plothandle=plot(geom.Z{k}, geom.R{k},'k-');
hold on
end
axis equal
[~,cont1]=contour(Z,R,R.*Aphi,15,'r:');
xlim([min(z) 0.04])
ylim([0.055 0.085])
%[~,cont2]=contour(Zphi,Rphi,Phi,20,'b');
rectangle('Position',[-0.011, 0.06375, 0.032+0.011, 0.081-0.06375],'EdgeColor','magenta','Linestyle','--')
t=linspace(0,-pi,500);
x=0.028/2*cos(t)+0.012;y=0.004*sin(t)+0.083;
phandle2=plot(x,y,'b-','displayname','Wall approximation 1');
x=0.028/2*cos(t)+0.012;y=0.012*sin(t)+0.091;
phandle3=plot(x,y,'g--','displayname','Wall approximation 2');
x=0.028/2*cos(t)+0.012;y=0.006*sin(t)+0.085;
phandle4=plot(x,y,'r--','displayname','Wall approximation 3');

legend([plothandle,cont1,phandle2,phandle3],{'Gun geometry', 'Magnetic field lines','Wall approximation'},'location','southwest')
f.PaperUnits='centimeters';
f.PaperSize=[12,8];
xlabel('z [m]')
ylabel('r [m]')
print(f,'phiBprofile','-dpdf','-fillpage')
savefig(f,'phiBprofile')
hold off




figure; contour(Z,R,Bz)
xlim([min(z) max(z)])
ylim([min(r) max(r)])

figure; contour(Z,R,Br)
xlim([min(z) max(z)])
ylim([min(r) max(r)])

function phi=Phivacuum(r,a,b,phia,phib)
    phi=((phib-phia)*log(r)+phia*log(b)-phib*log(a))/log(b/a);
end
