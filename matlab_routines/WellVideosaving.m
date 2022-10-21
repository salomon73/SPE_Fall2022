
fieldstep=1;

maxN=1.4e17;
Minwell=-6000;
%tend=min(8000,length(M.t2d));
tend=length(M.t2d);
plot3d=true;
logdensity=false;

fracn=0.1;

rAthetpos=177;

thefig=figure('Position',[0 0 1600 900]);
filename=M.name;
if logdensity
    filename=strcat(filename,'_log');
end
if plot3d
    filename=strcat(filename,'_well3D.avi');
else
    filename=strcat(filename,'_well2D.avi');
end
videowriterobj=VideoWriter([M.folder,'/',filename]);
videowriterobj.FrameRate=5;
open(videowriterobj);



if plot3d
    plotaxes(1)=subplot(2,1,1,'Parent',thefig);
    plotaxes(2)=subplot(2,1,2,'Parent',thefig);
    sf=surface(plotaxes(1),M.zgrid,M.rgrid,M.N(:,:,fieldstep),'edgecolor','none');
    
    hold(plotaxes(1),'on')
    contour(plotaxes(1),M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r--')
    if logdensity
            set(plotaxes(1),'zscale','log')
            set(plotaxes(1),'colorscale','log')
    end
    xlim(plotaxes(1),[M.zgrid(1) M.zgrid(end)])
    %ylim(plotaxes(1),[M.rgrid(1) M.rgrid(sum(M.nnr(1:2))+5)])
    xlabel(plotaxes(1),'z [m]')
    ylabel(plotaxes(1),'r [m]')
    title(plotaxes(1),'Density')
    c = colorbar(plotaxes(1));
    c.Label.String= 'n[m^{-3}]';
    %caxis(plotaxes(1),[-Inf maxN]);
    
    
    
    well=M.PotentialWell(fieldstep);
    sf=surface(plotaxes(2),M.zgrid,M.rgrid,well','edgecolor','none');
    hold(plotaxes(2),'on')
    dens=M.N(:,:,fieldstep);
    contour(plotaxes(2),M.zgrid,M.rgrid,dens-fracn*max(dens(:)),[0 0],'k--')
    contour(plotaxes(2),M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r--')
    
    xlim(plotaxes(2),[M.zgrid(1) M.zgrid(end)])
    %ylim(plotaxes(2),[M.rgrid(1) M.rgrid(sum(M.nnr(1:2))+5)])
    xlabel(plotaxes(2),'z [m]')
    ylabel(plotaxes(2),'r [m]')
    title(plotaxes(2),'Well')
    c = colorbar(plotaxes(2));
    c.Label.String= 'depth [V]';
    view(plotaxes(2),2)
    colormap(plotaxes(2),'jet');
    caxis(plotaxes(2),[Minwell 0]);
    zlim(plotaxes(2),[Minwell 0]);
    
    axis(plotaxes(1),'equal')
    axis(plotaxes(2),'equal')
    
    
    for i=fieldstep:80:tend
        sgtitle(thefig,sprintf('t= %1.2f [ns]',M.t2d(i)*1e9))
        dens=M.N(:,:,i);
        dens(M.geomweight(:,:,1)<0)=NaN;
        
        well=M.PotentialWell(i);
        well=well';
        
        plotaxes(1).Children(end).ZData=dens;
        plotaxes(1).Children(end).CData=dens;
        
        plotaxes(2).Children(end-1).ZData=dens-fracn*max(dens(:));
        plotaxes(2).Children(end).ZData=well;
        plotaxes(2).Children(end).CData=well;
        writeVideo(videowriterobj,getframe(thefig));
    end
else
    dens=M.N(:,:,fieldstep);
    model=M.potentialwellmodel(fieldstep);
    z=model.z;
    r=model.r;
    pot=model.pot;
    rathet=model.rathet;
    [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,1));
    dens=griddata(Zmesh,M.rAthet,dens,Zmesh,Rmesh);
    plot(M.zgrid,dens(rAthetpos,1:end));
    plotaxes(1)=gca;
    xlim(plotaxes(1),[M.zgrid(1) M.zgrid(end)])
    ylim(plotaxes(1),[0 maxN]);
    xlabel(plotaxes(1),'z [m]')
    ylabel(plotaxes(1),'n [m^{-3}]')
    grid(plotaxes(1), 'on')
    
    pot=griddata(z,rathet,pot,Zmesh,Rmesh);
    yyaxis right
    
    plot(M.zgrid(1:end),pot(rAthetpos,1:end));
    plotaxes(2)=gca;
    ylim(plotaxes(2),[Minwell 0]);
    ylabel(plotaxes(2),'depth [eV]')
    Curve2=plotaxes(2).Children;
    yyaxis left
    plotaxes(1)=gca;
    Curve1=plotaxes(1).Children;
    
    for i=fieldstep:80:tend
        sgtitle(thefig,sprintf('rA_\\theta=%1.3g [Tm^2] t= %1.2f [ns]',Rmesh(rAthetpos,1),M.t2d(i)*1e9))
        dens=M.N(:,:,i);
        model=M.potentialwellmodel(i);
        z=model.z;
        r=model.r;
        pot=model.pot;
        rathet=model.rathet;
        [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,1));
        dens=griddata(Zmesh,M.rAthet,dens,Zmesh,Rmesh);
        pot=griddata(z,rathet,pot,Zmesh,Rmesh);
        
        Curve1.YData=dens(rAthetpos,1:end);
        Curve2.YData=pot(rAthetpos,1:end);
        drawnow;
        writeVideo(videowriterobj,getframe(thefig));
    end
end
close(videowriterobj);