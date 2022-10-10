%geom=importgeometry('../geometry.data');

z=linspace(0.225,0.525,2500);
r=linspace(0.00,0.10,1000);
[Z,R] = meshgrid(z,r);


% 100 Amperes current in the coils
I = 100*ones(1,13); 
res=zeros([size(Z),3]);

for i=1:size(R,2)
res(:,i,:) = B_Ellip_10T_DNP({'aphi','bz','br'},'new',I,R(:,i),Z(:,i));
end
Aphi=res(:,:,1);
Bz=res(:,:,2);
Br=res(:,:,3);


% Plot the magnetic field lines
f=figure;
% for k=1:size(geom.Z)
% plothandle=plot(geom.Z{k}, geom.R{k},'k-');
% hold on
% end
axis equal
[~,cont1]=contour(Z,R,R.*Aphi,15,'r:');
%xlim([min(z) 0.04])
ylim([0.0 0.085])
%[~,cont2]=contour(Zphi,Rphi,Phi,20,'b');

f.PaperUnits='centimeters';
f.PaperSize=[12,8];
xlabel('z [m]')
ylabel('r [m]')
print(f,'phiBprofile_10T','-dpdf','-fillpage')
savefig(f,'phiBprofile_10T')
hold off




figure; contourf(Z,R,Bz)
hold on
contour(Z,R,R.*Aphi,15,'r:');

xlim([min(z) max(z)])
ylim([min(r) max(r)])


function phi=Phivacuum(r,a,b,phia,phib)
    phi=((phib-phia)*log(r)+phia*log(b)-phib*log(a))/log(b/a);
end
