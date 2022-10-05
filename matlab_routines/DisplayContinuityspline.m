function DisplayContinuityspline(M,Continuity,norm,log,contourscale, zlims)
%ForceBalance Show the radial force balance
%   Plot the three radial forces for the given time-step it or averaged
%   over the range of time steps defined in it

switch nargin
    case 2
        contourscale=0.2;
        norm=false;
        log=false;
        zlims=[-inf inf];
    case 3
        log=false;
        contourscale=0.2;
        zlims=[-inf inf];
    case 4
        contourscale=0.2;
        zlims=[-inf inf];
    case 5
        zlims=[-inf inf];
    case 6
        
    otherwise
        error("Invalid number of arguments")
end

contourcolor=[0 0 0];%[255 20 147]/255;
it=Continuity.it;

N=mean(Continuity.N,3);
maxdens=max(N(:));
densitycontour=contourc(M.zgrid,M.rgrid,mean(N,3),[contourscale contourscale]*maxdens);
rgridmax=sum(M.nnr(1:2));

S=1*mean(Continuity.S,3);


%%
f3=figure();
ax7=subplot(5,1,1);
zlim=[-1 1]*1e8;
plotvalue(ax7,mean(Continuity.dndt,3)./N,'\partialn/\partialt','\partialn/\partialt [m^{-3}s^{-1}]',true,zlim)


%zlim=[max([min(min(mean(sum(Continuity.ndivu(2:end,:,:),4),3))),min(min(mean(sum(Continuity.ugradn(2:end,:,:),4),3)))]) min([max(max(mean(sum(Continuity.ndivu(2:end,:,:),4),3))),max(max(mean(sum(Continuity.ugradn(2:end,:,:),4),3)))])];

ax8=subplot(5,1,2);
plotvalue(ax8,mean(sum(Continuity.ndivu,4),3)./N,'n\nabla(u)','n\nabla(u)[m^{-3}s^{-1}]',true,zlim)

ax9=subplot(5,1,3);
plotvalue(ax9,mean(sum(Continuity.ugradn,4),3)./N,'u\nabla(n)','u\nabla(n)[m^{-3}s^{-1}]',true,zlim)


ax10=subplot(5,1,4);
plotvalue(ax10,mean(Continuity.S,3)./N,'S_{io}','S_{io} [m^{-3}s^{-1}]',true,zlim)

ax11=subplot(5,1,5);
plotvalue(ax11,(mean(Continuity.dndt,3)+mean(sum(Continuity.ndivu,4)+sum(Continuity.ugradn,4),3)-S)./N,'Total','total[m^{-3}s^{-1}]',true,zlim)

linkaxes([ax7,ax8,ax9,ax10,ax11])
xlim([ax7,ax8,ax9,ax10,ax11],[1e-3 12e-3])
ylim([ax7,ax8,ax9,ax10,ax11],[0.078 0.0792])

sgtitle(sprintf('Continuity t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3}',M.t2d(min(it)),M.t2d(max(it)),double(maxdens)))
f3.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f3.PaperUnits='centimeters';
papsize=[20 19];
f3.PaperSize=papsize;
print(f3,sprintf('%scontinuity%d',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f3,sprintf('%scontinuity%d',name,floor(mean(it))))


%%
f3=figure();
ax7=subplot(5,1,1);
plotvalue(ax7,mean(Continuity.dndt,3)./N,'\partialn/\partialt','\partialn/\partialt [m^{-3}s^{-1}]',true,zlim)


%zlim=[max([min(min(mean(sum(Continuity.ndivu(2:end,:,:),4),3))),min(min(mean(sum(Continuity.ugradn(2:end,:,:),4),3)))]) min([max(max(mean(sum(Continuity.ndivu(2:end,:,:),4),3))),max(max(mean(sum(Continuity.ugradn(2:end,:,:),4),3)))])];

ax8=subplot(5,1,2);
plotvalue(ax8,mean(sum(Continuity.ndivu,4),3)./N,'n\nabla(u)','n\nabla(u)[m^{-3}s^{-1}]',true,zlim)

ax9=subplot(5,1,3);
plotvalue(ax9,mean(sum(Continuity.ugradn,4),3)./N,'u\nabla(n)','u\nabla(n)[m^{-3}s^{-1}]',true,zlim)


ax10=subplot(5,1,4);
plotvalue(ax10,mean(Continuity.S,3)./N,'S_{io}','S_{io} [m^{-3}s^{-1}]',true,zlim)

ax11=subplot(5,1,5);
plotvalue(ax11,(-mean(sum(Continuity.ndivu,4)+sum(Continuity.ugradn,4),3)+S)./N,'S-\nabla(nu)','S-\nabla(nu) [m^{-3}s^{-1}]',true,zlim)

linkaxes([ax7,ax8,ax9,ax10,ax11])
xlim([ax7,ax8,ax9,ax10,ax11],[1e-3 12e-3])
ylim([ax7,ax8,ax9,ax10,ax11],[0.078 0.0792])

sgtitle(sprintf('Continuity t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3}',M.t2d(min(it)),M.t2d(max(it)),double(maxdens)))
f3.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f3.PaperUnits='centimeters';
papsize=[20 19];
f3.PaperSize=papsize;
print(f3,sprintf('%scontinuity%d',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f3,sprintf('%scontinuity%d',name,floor(mean(it))))


%%
f3=figure;
ax6=subplot(4,1,1);
plotvalue(ax6,squeeze(mean(sum(Continuity.ndivu(:,:,:,:),4),3)),'ndivu er','radial ndivu [m^{-3}s^{-1}]',true,zlim)

ax7=subplot(4,1,2);
plotvalue(ax7,squeeze(mean(Continuity.ndivu(:,:,:,1),3)),'ndivu er1','radial 1 ndivu [m^{-3}s^{-1}]',true,zlim)

ax8=subplot(4,1,3);
plotvalue(ax8,squeeze(mean(Continuity.ndivu(:,:,:,2),3)),'ndivu er2','radial 2 ndivu [m^{-3}s^{-1}]',true,zlim)

ax9=subplot(4,1,4);
plotvalue(ax9,squeeze(mean(Continuity.ndivu(:,:,:,3),3)),'ndivu ez','axial ndivu [m^{-3}s^{-1}]',true,zlim)



linkaxes([ax6 ax7,ax8,ax9])
xlim([ax6 ax7,ax8,ax9],[1e-3 12e-3])
ylim([ax6 ax7,ax8,ax9],[0.0785 0.0792])

sgtitle(sprintf('Continuity ndivu t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3}',M.t2d(min(it)),M.t2d(max(it)),double(maxdens)))
f3.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f3.PaperUnits='centimeters';
papsize=[20 19];
f3.PaperSize=papsize;
print(f3,sprintf('%scontinuityndivu%d',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f3,sprintf('%scontinuityndivu%d',name,floor(mean(it))))
linkaxes([ax7,ax8,ax9])


%%
f3=figure;
ax7=subplot(2,1,1);
plotvalue(ax7,squeeze(mean(Continuity.ugradn(:,:,:,1),3)),'ugradn er','radial ugradn [m^{-3}s^{-1}]',true,zlim)

ax8=subplot(2,1,2);
plotvalue(ax8,squeeze(mean(Continuity.ugradn(:,:,:,2),3)),'ugradn ez','axial ugradn [m^{-3}s^{-1}]',true,zlim)


linkaxes([ax7,ax8])
xlim([ax7,ax8],[1e-3 12e-3])
ylim([ax7,ax8],[0.0785 0.0792])

sgtitle(sprintf('Continuity ugradn t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3}',M.t2d(min(it)),M.t2d(max(it)),double(maxdens)))
f3.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f3.PaperUnits='centimeters';
papsize=[20 19];
f3.PaperSize=papsize;
print(f3,sprintf('%scontinuityugradn%d',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f3,sprintf('%scontinuityugradn%d',name,floor(mean(it))))
linkaxes([ax7,ax8])


    function plotvalue(axis,matrix,titl, clab, rm1strow, zlims)
        if nargin < 5
            rm1strow=false;
        end
        if nargin>5
            zmin=zlims(1);
            zmax=zlims(2);
        else
            zmin=-inf;
            zmax=inf;
        end
        if rm1strow
            if log
                h=surface(axis,M.zgrid,M.rgrid(2:end),abs(matrix(2:end,:)));
            else
                h=surface(axis,M.zgrid,M.rgrid(2:end),matrix(2:end,:));
            end
        else
            if log
                h=surface(axis,M.zgrid,M.rgrid,abs(matrix));
            else
                h=surface(axis,M.zgrid,M.rgrid,matrix);
            end
        end
        %set(h,'edgecolor','none');
        hold on;
        if log
            zpos=interp2(M.zgrid, M.rgrid, abs(matrix), densitycontour(1,2:end), densitycontour(2,2:end));
        else
            zpos=interp2(M.zgrid, M.rgrid, matrix, densitycontour(1,2:end), densitycontour(2,2:end));
        end
        plot3(axis,densitycontour(1,2:end),densitycontour(2,2:end),zpos,'-.','linewidth',2,'color',contourcolor)
        xlabel(axis,'z [m]')
        ylabel(axis,'r [m]')
        xlim(axis,[M.zgrid(1) M.zgrid(end)])
        if M.conformgeom
            ylim(axis,[M.rgrid(1) M.rgrid(rgridmax)])
        end
        colormap(axis,'jet')
        c = colorbar(axis);
        c.Label.String=clab;
        caxis(axis,[zmin zmax])
        if norm
            caxis(axis,[Fminr Fmaxr]);
        end
        if log
            set(axis,'colorscale','log')
        end
        title(axis,titl)
        grid(axis, 'on')
        view(axis,2)
    end
end

