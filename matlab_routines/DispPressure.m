function DispPressure(M,it)

dens=mean(M.N(:,:,it),3);
maxdens=max(max(dens));

f=figure('title',M.name);
plotsubplot(subplot(2,3,1),squeeze(mean(M.Presstens(1,:,:,it),4)),'P_{rr} [Pa]');
plotsubplot(subplot(2,3,2),squeeze(mean(M.Presstens(2,:,:,it),4)),'P_{r\theta} [Pa]');
plotsubplot(subplot(2,3,3),squeeze(mean(M.Presstens(3,:,:,it),4)),'P_{rz} [Pa]');
plotsubplot(subplot(2,3,4),squeeze(mean(M.Presstens(4,:,:,it),4)),'P_{\theta\theta} [Pa]');
plotsubplot(subplot(2,3,5),squeeze(mean(M.Presstens(5,:,:,it),4)),'P_{z\theta} [Pa]');
plotsubplot(subplot(2,3,6),squeeze(mean(M.Presstens(6,:,:,it),4)),'P_{zz} [Pa]');

sgtitle(sprintf('Pressure t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(min(it)),M.t2d(max(it)),double(maxdens)))
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[16 14];
f.PaperSize=papsize;
print(f,sprintf('%sfluid_Pressure',name),'-dpdf','-fillpage')

function plotsubplot(ax,Val,clabel)
    h=surface(ax,M.zgrid,M.rgrid,Val);
    xlim(ax,[M.zgrid(1) M.zgrid(end)])
    rgridend=sum(M.nnr(1:2));
    ylim(ax,[M.rgrid(1) M.rgrid(rgridend)])
    xlabel(ax,'z [m]')
    ylabel(ax,'r [m]')
    c = colorbar(ax);
    c.Label.String= clabel;
    view(ax,2)
    set(h,'edgecolor','none');
end
end

