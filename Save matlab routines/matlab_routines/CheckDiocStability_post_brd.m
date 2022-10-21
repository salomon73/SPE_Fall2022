nbzid=15;
zindices=linspace(length(M.zgrid)/6,length(M.zgrid)*5/6,nbzid);
zconsid=2:2:nbzid-1;


dens=mean(M.N(:,:,fieldstart:fieldend),3);
wd=M.qe*dens./(2*M.eps_0*M.B');

wd=max(wd(:));
f=figure('Name', sprintf('%s dioc_stab',M.name));
fieldstart=max(1,length(M.t2d)-500);
fieldend=length(M.t2d);
ax1=subplot(2,1,1);
geomw=M.geomweight(:,:,1);
dens=mean(M.N(:,:,fieldstart:fieldend),3);
dens(geomw<=0)=0;
[C,h]=contourf(ax1,M.zgrid*1e3,M.rgrid*1e3,dens,'LineColor','none');
set(h,'LineColor','none')
hold on
contour(ax1,M.zgrid*1e3,M.rgrid*1e3,M.geomweight(:,:,1),[0 0],'w-.','linewidth',1.5,'Displayname','Boundaries');
if(M.conformgeom)
    ylim(ax1,[M.rgrid(1) M.rgrid(sum(M.nnr(1:2)))]*1e3)
else
    ylim(ax1,[M.rgrid(1) M.rgrid(end)]*1e3)
end
xlim(ax1,[M.zgrid(1) M.zgrid(end)]*1e3)
xlabel(ax1,'z [mm]')
ylabel(ax1,'r [mm]')
c = colorbar(ax1);
c.Label.String= 'n[m^{-3}]';
view(ax1,2)
grid on;
hold on;

zpos=zeros(size(zconsid));
wrslt=zeros(size(zconsid));
lmax=zeros(size(zconsid));
for zid=1:length(zconsid)%1:length(zindices)%;%floor(length(M.zgrid)/2);
    zindex=floor(zindices(zconsid(zid)));
    plot(ax1,M.zgrid(zindex)*[1 1]*1e3,[M.rgrid(1) M.rgrid(end)]*1e3,'r--')
    zpos(zid)=M.zgrid(zindex)*1e3;
    zfolder=sprintf('diocz_%03d',zindex);
    if ~isfolder(zfolder) || ~isfile(sprintf('%s/results.mat',zfolder))
        fprintf('Warning: result %s doesn''t exist\n',zfolder)
        continue
    end
    result=load(sprintf('%s/results.mat',zfolder));
    if ~isfield(result,'geomweight')
        result.geomweight=M.geomweight(:,:,1);
    end
    rindex=M.geomweight(:,zindex,1)>=0;
    
    for l=result.lsteps
        
        r=result.r(:,l);
        phin=result.phinum(:,:,l);
        w=result.w(:,l);
        [~,idwmax]=max(imag(w));
        if(imag(w(idwmax))>imag(wrslt(zid)))
            wrslt(zid)=w(idwmax);
            lmax(zid)=l;
        end
    end
end
ax2=subplot(2,1,2);
semilogy(zpos,real(wrslt)/wd,'bx','displayname','Real')
xlabel('z [mm]')
hold on
semilogy(zpos,imag(wrslt)/wd,'b+','displayname','Imag')
xlim([M.zgrid(1) M.zgrid(end)]*1e3)
ylabel('\omega/\omega_d')
ylim([1e-2 2])
yyaxis right
plot(zpos,lmax,'d','displayname','l')
ylabel('most unstable l')
grid on
title('Most unstable modes')
legend('location','east')


drawnow
posit1=get(ax1,'position');
posit2=get(ax2,'position');
posit2(3)=posit1(3);
set(ax2,'position',posit2);

f.PaperOrientation='landscape';
f.PaperUnits='centimeters';
papsize=[10 12];
f.PaperSize=papsize;
print(f,sprintf('%sdioc_stab',M.name),'-dpdf','-fillpage')
savefig(f,sprintf('%sdioc_stab',M.name))