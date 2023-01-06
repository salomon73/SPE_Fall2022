        
    fieldstep = length(ions.t2d);
    fixed = false;
    logdensity=false;
    showgrid=false;

        figure
        dens=ions.N(:,:,fieldstep);
        dens(ions.geomweight(:,:,1)<0)=NaN;
        %[~,sf]=contourf(ax1,ions.zgrid,ions.rgrid,dens,40,'edgecolor','none');
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
        %c.Limits=[0 max(ions.N(:))];
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
        
        Blines=ions.rAthet;
        levels=linspace(min(Blines(ions.geomweight(:,:,1)>0)),max(Blines(ions.geomweight(:,:,1)>0)),20);
        Blines(ions.geomweight(:,:,1)<0)=NaN;
        contour(ions.zgrid,ions.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
        set (gca, 'fontsize', 20)
        
        
        %%
        figure
            subplot(1,2,1)
                fieldstep = length(ions.t2d);
                fixed = false;
                logdensity=false;
                showgrid=false;
                dens=ions.N(:,:,fieldstep);
                dens(ions.geomweight(:,:,1)<0)=NaN;
                %[~,sf]=contourf(ax1,ions.zgrid,ions.rgrid,dens,40,'edgecolor','none');
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
                %c.Limits=[0 max(ions.N(:))];
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

                Blines=ions.rAthet;
                levels=linspace(min(Blines(ions.geomweight(:,:,1)>0)),max(Blines(ions.geomweight(:,:,1)>0)),20);
                Blines(ions.geomweight(:,:,1)<0)=NaN;
                contour(ions.zgrid,ions.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
                set (gca, 'fontsize', 20)
        subplot(1,2,2)
                fieldstep = length(ions_2.t2d);
                fixed = false;
                logdensity=false;
                showgrid=false;
                dens=ions_2.N(:,:,fieldstep);
                dens(ions_2.geomweight(:,:,1)<0)=NaN;
                %[~,sf]=contourf(ax1,ions.zgrid,ions.rgrid,dens,40,'edgecolor','none');
                if(showgrid)
                [sf]=surface(ions_2.zgrid,ions_2.rgrid,dens);
                else
                [~,sf]=contourf(ions_2.zgrid,ions_2.rgrid,dens,40,'edgecolor','none');
                end
                xlim([ions_2.zgrid(1) ions_2.zgrid(end)])
                ylim([ions_2.rgrid(1) ions_2.rgrid(end)])
                xlabel('z [m]')
                ylabel('r [m]')
                title('Density')
                c = colorbar;
                c.Label.String= 'n[m^{-3}]';
                %c.Limits=[0 max(ions.N(:))];
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

                contour(ions_2.zgrid,ions_2.rgrid,ions_2.geomweight(:,:,1),[0 0],'r-','linewidth',1.5);

                Blines=ions.rAthet;
                levels=linspace(min(Blines(ions_2.geomweight(:,:,1)>0)),max(Blines(ions_2.geomweight(:,:,1)>0)),20);
                Blines(ions_2.geomweight(:,:,1)<0)=NaN;
                contour(ions_2.zgrid,ions_2.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
                set (gca, 'fontsize', 20)
