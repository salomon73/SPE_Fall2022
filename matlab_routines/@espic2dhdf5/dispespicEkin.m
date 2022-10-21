function f=dispespicEkin(M)
fieldstep=1;


f=uifigure('Name','Fluid Energy data');

mf=uipanel(f,'Position',[5 50 f.Position(3)-10 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);

sgtitle(mf,sprintf('t=%0.5e s',M.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.6*m.Position(3) 3]);
sld.Value=fieldstep;
sld.Limits=[1 length(M.t2d)];

edt = uieditfield(m,'numeric','Limits',[1 length(M.t2d)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';


Printbt=uibutton(m,'Position',[edt.Position(1)+edt.Position(3)+10 5 40 20],'Text', 'Save');
%Playbt=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play/Pause');


sld.ValueChangingFcn={@updatefigdata,edt,mf};
edt.ValueChangedFcn={@updatefigdata,sld,mf};
Printbt.ButtonPushedFcn={@plotGridButtonPushed};
[R,Z]=meshgrid(M.rgrid,M.zgrid);
Rinv=1./R;
Rinv(:,1)=0;

PlotEspic2ddata(mf,M,fieldstep);

    function PlotEspic2ddata(fig,M,fieldstep)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        sgtitle(fig,sprintf('t=%0.5e s',M.t2d(fieldstep)))
        
        
        ax1=subplot(2,2,1,'Parent',fig);
        
        
        surface(ax1,M.zgrid,M.rgrid,M.N(:,:,fieldstep),'edgecolor','none');
        hold(ax1,'on')
        
        border=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);
        zpos=interp2(M.zgrid, M.rgrid, M.N(:,:,fieldstep), border(1,2:end), border(2,2:end));
        plot3(ax1,border(1,2:end),border(2,2:end),zpos,'r-','linewidth',1.5,'Displayname','Boundaries')
        
        bfields=contourc(M.zgrid,M.rgrid,M.rAthet,15);
        zpos=interp2(M.zgrid, M.rgrid, M.N(:,:,fieldstep), bfields(1,2:end), bfields(2,2:end));
        plot3(ax1,bfields(1,2:end),bfields(2,2:end),zpos,'m-.','linewidth',1.5,'Displayname','Boundaries')
        
        
        xlim(ax1,[M.zgrid(1) M.zgrid(end)])
        ylim(ax1,[M.rgrid(1) M.rgrid(end)])
        xlabel(ax1,'z [m]')
        ylabel(ax1,'r [m]')
        title(ax1,'Position')
        c = colorbar(ax1);
        c.Label.String= 'n[m^{-3}]';
        %c.Limits=[0 max(M.N(:))];
        %caxis(ax1,[0 max(M.N(:))]);
        view(ax1,2)
        
        %set(ax1,'colorscale','log')
        
        
        UR=M.fluidUR(:,:,fieldstep);
        UZ=M.fluidUZ(:,:,fieldstep);
        
        ax2=subplot(2,2,2,'Parent',fig);
        
        surface(ax2,M.zgrid,M.rgrid,squeeze(M.fluidEkin(1,:,:,fieldstep))/M.qe,'edgecolor','none');
        %plot(ax2,M.zgrid,data.pot(:,5))
        xlabel(ax2,'z [m]')
        ylabel(ax2,'r [m]')
        colormap(ax2,'jet')
        c = colorbar(ax2);
        c.Label.String= 'E_r [eV]';
        hold(ax2,'on')
        border=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);
        zpos=interp2(M.zgrid, M.rgrid, squeeze(M.fluidEkin(1,:,:,fieldstep))/M.qe, border(1,2:end), border(2,2:end));
        plot3(ax2,border(1,2:end),border(2,2:end),zpos,'r-','linewidth',1.5,'Displayname','Boundaries')
        %c.Limits=[min(M.fluidUR(:,:,:)) max(M.fluidUR(:,:,:))];
        
        grid(ax2, 'on')
        %caxis(ax2,[-4e5 1e6])
        
        view(ax2,2)
        
        ax3=subplot(2,2,3,'Parent',fig);
        
        surface(ax3,M.zgrid,M.rgrid,squeeze(M.fluidEkin(2,:,:,fieldstep))/M.qe,'edgecolor','none')
        
        xlabel(ax3,'z [m]')
        ylabel(ax3,'r [m]')
        colormap(ax3,'jet')
        c = colorbar(ax3);
        c.Label.String= 'E_\theta [eV]';
        hold(ax3,'on')
        border=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);
        zpos=interp2(M.zgrid, M.rgrid, squeeze(M.fluidEkin(2,:,:,fieldstep))/M.qe, border(1,2:end), border(2,2:end));
        plot3(ax3,border(1,2:end),border(2,2:end),zpos,'r-','linewidth',1.5,'Displayname','Boundaries')

        
        grid(ax3, 'on')
        view(ax3,2)
        
        ax4=subplot(2,2,4,'Parent',fig);
        
        surface(ax4,M.zgrid,M.rgrid,squeeze(M.fluidEkin(3,:,:,fieldstep))/M.qe,'edgecolor','none')
        xlabel(ax4,'z [m]')
        ylabel(ax4,'r [m]')
        colormap(ax4,'jet')
        c = colorbar(ax4);
        c.Label.String= 'E_z [eV]';
        titl='';
        labl='';
        hold(ax4,'on')
        border=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);
        zpos=interp2(M.zgrid, M.rgrid, squeeze(M.fluidEkin(3,:,:,fieldstep))/M.qe, border(1,2:end), border(2,2:end));
        plot3(ax4,border(1,2:end),border(2,2:end),zpos,'r-','linewidth',1.5,'Displayname','Boundaries')

        
        grid(ax4, 'on')
        view(ax4,2)
        
        linkaxes([ax1 ax2 ax3 ax4],'xy')
        
    end

    function plotGridButtonPushed(btn,ax)
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        f=figure();
        PlotEspic2dgriddata(f,M,sld.Value);
        f.PaperOrientation='landscape';
        [~, name, ~] = fileparts(M.file);
        print(f,sprintf('%sfluid%d',name,sld.Value),'-dpdf','-fillpage')
    end

    function updatefigdata(control, event, Othercontrol, fig)
        
        
        if strcmp(event.EventName,'ValueChanged')
            fieldstep=floor(control.Value);
            control.Value=fieldstep;
        else
            fieldstep=floor(event.Value);
        end
        Othercontrol.Value=fieldstep;
        
        sgtitle(fig,sprintf('t=%0.5e s',double((fieldstep-1)*M.it1)*M.dt))
        
        %% update Position histogram
        ax1=fig.Children(end);
        
        
        dens=M.N(:,:,fieldstep);
        dens(M.geomweight(:,:,1)<0)=0;

        zpos=interp2(M.zgrid, M.rgrid, dens, ax1.Children(end-1).XData, ax1.Children(end-1).YData);
        ax1.Children(end-1).ZData=zpos;
        zpos=interp2(M.zgrid, M.rgrid, dens, ax1.Children(end-2).XData, ax1.Children(end-2).YData);
        ax1.Children(end-2).ZData=zpos;
        
        ax1.Children(end).ZData=dens;
        ax1.Children(end).CData=dens;
        
        
        ER=squeeze(M.fluidEkin(1,:,:,fieldstep))/M.qe;
        EZ=squeeze(M.fluidEkin(3,:,:,fieldstep))/M.qe;
        
        view(ax1,2)
        %% update Radial velocity
        ax1=fig.Children(end-2);
        zpos=interp2(M.zgrid, M.rgrid, ER, ax1.Children(end-1).XData, ax1.Children(end-1).YData);
        ax1.Children(end-1).ZData=zpos;
        ER(ER<=0)=NaN;
        fig.Children(end-2).Children(end).CData=ER;
        fig.Children(end-2).Children(end).ZData=ER;
        caxis(ax1,[0 50]);
        
        
        %% update Azimuthal velocity
        ax1=fig.Children(end-4);
        Ethet=squeeze(M.fluidEkin(2,:,:,fieldstep))/M.qe;
        zpos=interp2(M.zgrid, M.rgrid, Ethet, ax1.Children(end-1).XData, ax1.Children(end-1).YData);
        ax1.Children(end-1).ZData=zpos;
        Ethet(Ethet<=0)=NaN;     
        fig.Children(end-4).Children(end).CData=Ethet;
        fig.Children(end-4).Children(end).ZData=Ethet;
        caxis(ax1,[0 50]);
        
        %% update Axial velocity
        ax1=fig.Children(end-6);
        zpos=interp2(M.zgrid, M.rgrid, EZ, ax1.Children(end-1).XData, ax1.Children(end-1).YData);
        ax1.Children(end-1).ZData=zpos;
        %drawnow limitrate
        EZ(EZ<=0)=NaN;
        fig.Children(end-6).Children(end).CData=EZ;
        fig.Children(end-6).Children(end).ZData=EZ;
        caxis(ax1,[0 50]);
        
        
    end

end


