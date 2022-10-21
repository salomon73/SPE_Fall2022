function dispespicFields(obj,logdensity,showgrid,fixed,parper)
%dispespicFields Allows to display the time evolution of the density, electric potential and electric fields
%   M is of class espic2dhdf5 and contains the simulation results
fieldstep=1;
if nargin <2
    logdensity=false;
    showgrid=false;
    fixed=false;
end
if nargin <3
    showgrid=false;
    fixed=false;
end
if nargin <4
    fixed=false;
end
if nargin <5
    parper=false;
end
fixed=fi(fixed);

f=uifigure('Name',sprintf('Grid data %s',obj.name));

mf=uipanel(f,'Position',[5 50 f.Position(3)-10 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);

sgtitle(mf,sprintf('step=%d t=%0.5e s',fieldstep*obj.it1,obj.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.6*m.Position(3) 3]);
sld.Value=fieldstep;
sld.Limits=[1 size(obj.t2d,1)];

edt = uieditfield(m,'numeric','Limits',[1 size(obj.t2d,1)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';

MaxN=0;
Printbt=uibutton(m,'Position',[edt.Position(1)+edt.Position(3)+10 5 40 20],'Text', 'Save');
Play=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play');
Pause=uibutton(m,'Position',[Play.Position(1)+Play.Position(3)+10 5 40 20],'Text', 'Pause');
%Playbt=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play/Pause');

stop=false;

sld.ValueChangingFcn={@updatefigdata,edt,mf};
edt.ValueChangedFcn={@updatefigdata,sld,mf};
Printbt.ButtonPushedFcn={@plotGridButtonPushed};
Play.ButtonPushedFcn={@plotPlayButtonPushed};

Pause.ButtonPushedFcn={@PauseButtonPushed};

set(f,'KeyPressFcn',{ @onKeyDown,sld,edt,mf})

PlotEspic2dgriddata(mf,obj,fieldstep);

    function plotPlayButtonPushed(btn,ax)
        stop=false;
        i=sld.Value;
        while ~stop
            edt.Value=i;
            sld.Value=i;
            updatesubplotsdata(i,mf);
            pause(0.01)
            i=sld.Value;
            i=i+10;
            if(i>sld.Limits(2))
                stop=true;
            end
        end
    end
    function PauseButtonPushed(btn,ax)
        stop = true;
    end

    function onKeyDown(src,event,slider,editfield, fig)
        direction=0;
        if strcmp(event.Key,'leftarrow')
            direction=-1;
        elseif strcmp(event.Key,'rightarrow')
            direction=+1;
        elseif strcmp(event.Key,'uparrow')
            direction=+10;
        elseif strcmp(event.Key,'downarrow')
            direction=-10;
        end
    
        if(direction~=0)
            currval=slider.Value;
            slider.Value=max(slider.Limits(1),min(currval+direction,slider.Limits(2)));
            updatefigdata(slider, event, editfield ,fig)
        end
    end

    function PlotEspic2dgriddata(fig,M,fieldstep)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*M.it1,M.t2d(fieldstep)))
        
        ax1=subplot(2,2,1,'Parent',fig);
        %%
        dens=obj.N(:,:,fieldstep);
        dens(obj.geomweight(:,:,1)<0)=NaN;
        %[~,sf]=contourf(ax1,M.zgrid,M.rgrid,dens,40,'edgecolor','none');
        if(showgrid)
            [sf]=surface(ax1,M.zgrid,M.rgrid,dens);
        else
            [~,sf]=contourf(ax1,M.zgrid,M.rgrid,dens,40,'edgecolor','none');
        end
        xlim(ax1,[M.zgrid(1) M.zgrid(end)])
        ylim(ax1,[M.rgrid(1) M.rgrid(end)])
        xlabel(ax1,'z [m]')
        ylabel(ax1,'r [m]')
        title(ax1,'Density')
        c = colorbar(ax1);
        c.Label.String= 'n[m^{-3}]';
        %c.Limits=[0 max(M.N(:))];
        if(isboolean(fixed) && fixed)
            climits=caxis(ax1);
            MaxN=climits(2);
        elseif(~isboolean(fixed))
            MaxN=fixed.data;
            caxis(ax1,[-Inf MaxN]);
        end
        view(ax1,2)
        hotmap=flipud(hot);
        
        colormap(ax1,hotmap);
        hold(ax1, 'on')
        %border=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);
        %zpos=interp2(M.zgrid, M.rgrid, M.N(:,:,fieldstep), border(1,2:end), border(2,2:end));
        %plot3(ax1,border(1,2:end),border(2,2:end),zpos+1e-6,'r-','linewidth',1.5,'Displayname','Boundaries')
        contour(ax1,M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r-','linewidth',1.5);
        
        Blines=M.rAthet;
        levels=linspace(min(Blines(M.geomweight(:,:,1)>0)),max(Blines(M.geomweight(:,:,1)>0)),20);
        Blines(M.geomweight(:,:,1)<0)=NaN;
%         bfields=contourc(M.zgrid,M.rgrid,Blines,real(levels));
%         zpos=interp2(M.zgrid, M.rgrid, M.N(:,:,fieldstep), bfields(1,2:end), bfields(2,2:end));
%         plot3(ax1,bfields(1,2:end),bfields(2,2:end),zpos+1e-6,'m-.','linewidth',1.5,'Displayname','Boundaries')
        
        contour(ax1,M.zgrid,M.rgrid,Blines,real(levels),'m-.','linewidth',1.5);
        
        if logdensity
            %set(ax1,'zscale','log')
            set(ax1,'colorscale','log')
        end
        
        ax2=subplot(2,2,2,'Parent',fig);
        pot=M.pot(:,:,fieldstep);
        pot(M.geomweight(:,:,1)<0)=NaN;
        contourf(ax2,M.zgrid,M.rgrid,pot,20);
        %plot(ax2,M.zgrid,data.pot(:,5))
        %%
        xlabel(ax2,'z [m]')
        ylabel(ax2,'r [m]')
        colormap(ax2,'jet')
        c = colorbar(ax2);
        c.Label.String= '\Phi [V]';
        title(ax2,'Es potential')
        grid(ax2, 'on')
        hold(ax2, 'on')
        contour(ax2,M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r-','linewidth',2.5,'Displayname','Boundaries');
        Blines=M.rAthet;
        levels=linspace(min(Blines(M.geomweight(:,:,1)>0)),max(Blines(M.geomweight(:,:,1)>0)),20);
        Blines(M.geomweight(:,:,1)<0)=NaN;
        contour(ax2,M.zgrid,M.rgrid,Blines,real(levels),'m-.','linewidth',1.5,'Displayname','Magnetic field lines');
        
        %%
        ax3=subplot(2,2,3,'Parent',fig);
        %%
        if parper
            Ez=M.Epar(fieldstep);
        else
            Ez=M.Ez(:,:,fieldstep);
        end
        Ez(M.geomweight(:,:,1)<0)=NaN;
        contourf(ax3,M.zgrid,M.rgrid,Ez,60)
        hold(ax3, 'on')
        contour(ax3,M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r-','linewidth',1.5,'Displayname','Boundaries');
        xlabel(ax3,'z [m]')
        ylabel(ax3,'r [m]')
        c = colorbar(ax3);
        c.Label.String= 'E_{z} [V/m]';
        if parper
            title(ax3,'Parallel Electric field')
        else
            title(ax3,'Axial Electric field')
        end
        grid(ax3, 'on')
        
        ax4=subplot(2,2,4,'Parent',fig);
        %%
        if parper
            Er=M.Eperp(fieldstep);
        else
            Er=M.Er(:,:,fieldstep);
        end
        Er(M.geomweight(:,:,1)<0)=NaN;
        contourf(ax4,M.zgrid,M.rgrid,Er,60)
        hold(ax4, 'on')
        contour(ax4,M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r-','linewidth',1.5,'Displayname','Boundaries');
        xlabel(ax4,'z [m]')
        ylabel(ax4,'r [m]')
        c = colorbar(ax4);
        c.Label.String= 'E_{r} [V/m]';
        if parper
            title(ax4,'Perpendicular Electric field')
        else
            title(ax4,'Radial Electric field')
        end
        grid(ax4, 'on')
        
        linkaxes([ax1,ax2,ax3,ax4]);
    end

    function plotGridButtonPushed(btn,ax)
        f=figure();
        PlotEspic2dgriddata(f,obj,floor(sld.Value));
        axold=mf.Children(end-2);
        xlimits=xlim(axold);
        ylimits=ylim(axold);
        
        ax1=subplot(2,2,1,'Parent',f);
        xlim(ax1,xlimits);
        ylim(ax1,ylimits);
        f.PaperOrientation='landscape';
        f.PaperUnits='centimeters';
        f.PaperSize=[18,14];
        [~, name, ~] = fileparts(obj.file);
        savefig(f,sprintf('%sGridFields%d',name,floor(sld.Value)))
        print(f,sprintf('%sGridFields%d',name,floor(sld.Value)),'-dpdf','-fillpage')
        set(f, 'Color', 'w');
        export_fig(f,sprintf('%sGridFields%d',name,floor(sld.Value)),'-eps')
    end

    function updatefigdata(control, event, Othercontrol, fig)
        
        
        if strcmp(event.EventName,'ValueChanged')
            fieldstep=floor(control.Value);
            control.Value=fieldstep;
        elseif strcmp(event.EventName,'KeyPress')
            fieldstep=floor(control.Value);
            control.Value=fieldstep;
        else
            fieldstep=floor(event.Value);
        end
        Othercontrol.Value=fieldstep;
        updatesubplotsdata(fieldstep, fig);
    end

    function updatesubplotsdata(fieldstep, fig)
        
        sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*obj.it1,obj.t2d(floor(fieldstep))))
        
        %% update Position histogram
        ax1=fig.Children(end);
        dens=obj.N(:,:,fieldstep);
        dens(obj.geomweight(:,:,1)<0)=NaN;
        if logdensity
            dens(dens<=0)=NaN;
        end
        lvls=linspace(0,max(dens(:)),40);
        xlimits=xlim(ax1);
        ylimits=ylim(ax1);
        
        %zpos=interp2(obj.zgrid, obj.rgrid, dens, ax1.Children(end-1).XData, ax1.Children(end-1).YData);
        %ax1.Children(end-1).ZData=zpos+1e-6;
        %zpos=interp2(obj.zgrid, obj.rgrid, dens, ax1.Children(end-2).XData, ax1.Children(end-2).YData);
        %ax1.Children(end-2).ZData=zpos+1e-6;
        
        ax1.Children(end).ZData=dens;
        ax1.Children(end).LevelList=lvls;
        if(showgrid)
            ax1.Children(end).CData=dens;
        end
        if(isboolean(fixed) && fixed)
            climits=caxis(ax1);
            MaxN=climits(2);
            nmax=max(dens(:));
            if(nmax>MaxN)
                MaxN=nmax;
            end
            caxis(ax1,[0 MaxN]);
        elseif(~isboolean(fixed))
            MaxN=fixed.data;
            caxis(ax1,[-Inf MaxN]);
        end
        caxis(ax1,'auto')
        
        %     view(ax1,2)
        %% update ES Potential
        
        ax2=subplot(2,2,2,'Parent',fig);
        %hold(ax2, 'off')
        pot=obj.pot(:,:,fieldstep);
        pot(obj.geomweight(:,:,1)<0)=NaN;
        contourf(ax2,obj.zgrid,obj.rgrid,pot,20);
        %ax2.Children(end).ZData=pot;
        caxis(ax2,[-inf inf])
        
%         xlabel(ax2,'z [m]')
%         ylabel(ax2,'r [m]')
%         colormap(ax2,'jet')
%         c = colorbar(ax2);
%         c.Label.String= '\Phi [V]';
%         title(ax2,'Es potential')
%         %%
%         hold(ax2, 'on')
%         contour(ax2,obj.zgrid,obj.rgrid,obj.geomweight(:,:,1),[0 0],'r-','linewidth',2.5,'Displayname','Boundaries');
%         Blines=obj.rAthet;
%         levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
%         Blines(obj.geomweight(:,:,1)<0)=NaN;
%         contour(ax2,obj.zgrid,obj.rgrid,Blines,real(levels),'m-.','linewidth',1.5,'Displayname','Magnetic field lines');
        
%         pot=obj.pot(:,:,fieldstep);
%         pot(obj.geomweight(:,:,1)<0)=0;
%         contourf(fig.Children(end-2),obj.zgrid,obj.rgrid,pot,20);


        %% update Radial electric field
        if parper
            Ez=obj.Epar(fieldstep);
        else
            Ez=obj.Ez(:,:,fieldstep);
        end
        Ez(obj.geomweight(:,:,1)<0)=NaN;
        ax3=subplot(2,2,3,'Parent',fig);
        ax3.Children(end).ZData=Ez;
        caxis(ax3,[-inf inf])
        %% update Axial electric field
        if parper
            Er=obj.Eperp(fieldstep);
        else
            Er=obj.Er(:,:,fieldstep);
        end
        Er(obj.geomweight(:,:,1)<0)=NaN;
        ax4=subplot(2,2,4,'Parent',fig);
        ax4.Children(end).ZData=Er;
        caxis(ax4,[-inf inf])
        
        xlim(ax4,xlimits);
        ylim(ax4,ylimits);
        %drawnow limitrate
        
    end

end


