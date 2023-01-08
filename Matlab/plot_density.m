        
    fieldstep = length(ions_less.t2d);
    fixed = false;
    logdensity=false;
    showgrid=false;

        figure
        dens=ions_less.N(:,:,fieldstep);
        dens(ions_less.geomweight(:,:,1)<0)=NaN;
        %[~,sf]=contourf(ax1,ions_less.zgrid,ions_less.rgrid,dens,40,'edgecolor','none');
        if(showgrid)
            [sf]=surface(ions_less.zgrid,ions_less.rgrid,dens);
        else
            [~,sf]=contourf(1000*ions_less.zgrid,1000*ions_less.rgrid,dens,40,'edgecolor','none');
        end
        xlim(1000*[ions_less.zgrid(1) ions_less.zgrid(end)])
        ylim(1000*[ions_less.rgrid(1) ions_less.rgrid(end)])
        xlabel('z [mm]')
        ylabel('r [mm]')
        title('Density')
        c = colorbar;
        c.Label.String= 'n[m^{-3}]';
        %c.Limits=[0 max(ions_less.N(:))];
        if(fixed)
            climits=caxis;
            MaxN=climits(2);
%         elseif(~isboolean(fixed))
%             MaxN=fixed.data;
%             caxis([-Inf MaxN]);
        end
        view(2)
        hotmap=flipud(hot);
        
        colormap(hotmap);
        hold on

        contour(1000*ions_less.zgrid,1000*ions_less.rgrid,ions_less.geomweight(:,:,1),[0 0],'r-','linewidth',1.5);
        
        Blines=ions_less.rAthet;
        levels=linspace(min(Blines(ions_less.geomweight(:,:,1)>0)),max(Blines(ions_less.geomweight(:,:,1)>0)),20);
        Blines(ions_less.geomweight(:,:,1)<0)=NaN;
        contour(1000*ions_less.zgrid,1000*ions_less.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
        set (gca, 'fontsize', 20)
        
        
        %%
        figure
            subplot(1,2,1)
                fieldstep = length(ions_less.t2d);
                fixed = false;
                logdensity=false;
                showgrid=false;
                dens=ions_less.N(:,:,fieldstep);
                dens(ions_less.geomweight(:,:,1)<0)=NaN;
                %[~,sf]=contourf(ax1,ions_less.zgrid,ions_less.rgrid,dens,40,'edgecolor','none');
                if(showgrid)
                [sf]=surface(ions_less.zgrid,ions_less.rgrid,dens);
                else
                [~,sf]=contourf(ions_less.zgrid,ions_less.rgrid,dens,40,'edgecolor','none');
                end
                xlim([ions_less.zgrid(1) ions_less.zgrid(end)])
                ylim([ions_less.rgrid(1) ions_less.rgrid(end)])
                xlabel('z [m]')
                ylabel('r [m]')
                title('Density')
                c = colorbar;
                c.Label.String= 'n[m^{-3}]';
                %c.Limits=[0 max(ions_less.N(:))];
                if(fixed)
                climits=caxis;
                MaxN=climits(2);
                %         elseif(~isboolean(fixed))
                %             MaxN=fixed.data;
                %             caxis([-Inf MaxN]);
                end
                view(2)
                hotmap=flipud(hot);

                colormap(hotmap);
                hold on

                contour(ions_less.zgrid,ions_less.rgrid,ions_less.geomweight(:,:,1),[0 0],'r-','linewidth',1.5);

                Blines=ions_less.rAthet;
                levels=linspace(min(Blines(ions_less.geomweight(:,:,1)>0)),max(Blines(ions_less.geomweight(:,:,1)>0)),20);
                Blines(ions_less.geomweight(:,:,1)<0)=NaN;
                contour(ions_less.zgrid,ions_less.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
                set (gca, 'fontsize', 20)
        subplot(1,2,2)
                fieldstep = length(ions.t2d);
                fixed = false;
                logdensity=false;
                showgrid=false;
                dens=ions.N(:,:,fieldstep);
                dens(ions.geomweight(:,:,1)<0)=NaN;
                %[~,sf]=contourf(ax1,ions_less.zgrid,ions_less.rgrid,dens,40,'edgecolor','none');
                if(showgrid)
                [sf]=surface(ions.zgrid,ions.rgrid,dens);
                else
                [~,sf]=contourf(ions.zgrid,ions.rgrid,dens,40,'edgecolor','none');
                end
                xlim([ions.zgrid(1) ions.zgrid(end)])
                ylim([ions.rgrid(1) ions.rgrid(end)])
                xlabel('z [m]')
                ylabel('r [m]')
                title('Density')
                c = colorbar;
                c.Label.String= 'n[m^{-3}]';
                %c.Limits=[0 max(ions_less.N(:))];
                if(fixed)
                climits=caxis;
                MaxN=climits(2);
                %         elseif(~isboolean(fixed))
                %             MaxN=fixed.data;
                %             caxis([-Inf MaxN]);
                end
                view(2)
                hotmap=flipud(hot);

                colormap(hotmap);
                hold on

                contour(ions.zgrid,ions.rgrid,ions.geomweight(:,:,1),[0 0],'r-','linewidth',1.5);

                Blines=ions_less.rAthet;
                levels=linspace(min(Blines(ions.geomweight(:,:,1)>0)),max(Blines(ions.geomweight(:,:,1)>0)),20);
                Blines(ions.geomweight(:,:,1)<0)=NaN;
                contour(ions.zgrid,ions.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
                set (gca, 'fontsize', 20)
