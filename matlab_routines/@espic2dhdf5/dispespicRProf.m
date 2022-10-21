function dispespicRProf(M,logdensity,showgrid,fixed,parper)
%dispespicFields Allows to display the time evolution of the density, electric potential and electric fields
%   M is of class espic2dhdf5 and contains the simulation results
fieldstep=1;
zpos=1;
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

f=uifigure('Name',sprintf('Grid data %s',M.name));

mf=uipanel(f,'Position',[5 50 f.Position(3)-30 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);
rpanel=uipanel(f,'Position',[f.Position(3)-28 50 25 f.Position(4)-55]);

sgtitle(mf,sprintf('step=%d t=%0.5e s',fieldstep*M.it1,M.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.6*m.Position(3) 3]);
sld.Value=fieldstep;
sld.Limits=[1 size(M.t2d,1)];
sld.Tag='timeslider';

sldr = uislider(rpanel,'Orientation','vertical','Position',[5 35 40 rpanel.Position(4)-40]);
sldr.Value=zpos;
sldr.Limits=[1 length(M.zgrid)];
sldr.Tag='axialslider';

edr = uieditfield(rpanel,'numeric','Limits',[1 length(M.zgrid)],'Value',1);
edr.Position=[sldr.Position(1) sldr.Position(2)-30 40 20];
edr.RoundFractionalValues='on';
edr.Tag='axialfield';

edt = uieditfield(m,'numeric','Limits',[1 size(M.t2d,1)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';
edt.Tag='timefield';

MaxN=0;
Printbt=uibutton(m,'Position',[edt.Position(1)+edt.Position(3)+10 5 40 20],'Text', 'Save');
Play=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play');
Pause=uibutton(m,'Position',[Play.Position(1)+Play.Position(3)+10 5 40 20],'Text', 'Pause');
%Playbt=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play/Pause');

stop=false;

sld.ValueChangingFcn={@updatefigdata,edt,mf};
edt.ValueChangedFcn={@updatefigdata,sld,mf};

sldr.ValueChangingFcn={@updatefigdata,edr,mf};
edr.ValueChangedFcn={@updatefigdata,sldr,mf};

Printbt.ButtonPushedFcn={@plotGridButtonPushed};
Play.ButtonPushedFcn={@plotPlayButtonPushed};

Pause.ButtonPushedFcn={@PauseButtonPushed};

set(f,'KeyPressFcn',{ @onKeyDown,sld,edt,mf})

PlotEspic2dgriddata(mf,M,fieldstep,zpos);

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

    function PlotEspic2dgriddata(fig,M,fieldstep,zpos)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        sgtitle(fig,sprintf('step=%d t=%0.5e s z=%0.3e mm',(fieldstep-1)*M.it1,M.t2d(fieldstep),M.zgrid(zpos)*1e3))
        
        geomw=M.geomweight(:,zpos,1)>=0;
        ax1=subplot(2,2,1,'Parent',fig);
        p=plot(ax1,M.rgrid*1e3,M.N(:,zpos,fieldstep),'linewidth',1.5);
        xlim(ax1,[M.rgrid(1) M.rgrid(end)]*1e3)
        xlabel(ax1,'r [mm]')
        title(ax1,'Density')
        ylabel(ax1,'n[m^{-3}]');
        %c.Limits=[0 max(M.N(:))];
        
        hold(ax1, 'on')
        [~,id1]=min(abs(M.geomweight(1:10,zpos,1)));
        [~,id2]=min(abs(M.geomweight(11:end,zpos,1)));
        id2=id2+10;
        ylimits=ylim;
        plot(ax1,M.rgrid(id1)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        plot(ax1,M.rgrid(id2)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        
        yyaxis(ax1,'right')
        hold(ax1, 'on')
        Er=M.Er(:,zpos,fieldstep).*geomw;
        Ez=M.Ez(:,zpos,fieldstep).*geomw;
        p1=plot(ax1,M.rgrid*1e3,Er,'linewidth',1.5);
        p2=plot(ax1,M.rgrid*1e3,Ez,'linewidth',1.5);
        ylabel(ax1,'E [V/m]')
        if max(abs([Er(:); Ez(:)]))>0 
            ylim(ax1,[ -max(abs([Er(:); Ez(:)])) max(abs([Er(:); Ez(:)]))])
        end
        legend(ax1,[p p1 p2],{'n','Er','Ez'},'location','northwest')

       
        
        ax2=subplot(2,2,2,'Parent',fig);
        ur=M.fluidUR(:,zpos,fieldstep);
        plot(ax2,M.rgrid*1e3,ur,'linewidth',1.5);
        xlim(ax2,[M.rgrid(1) M.rgrid(end)]*1e3)
        xlabel(ax2,'r [mm]')
        title(ax2,'radial velocity')
        ylabel(ax2,'v_r [m/s]');
        %c.Limits=[0 max(M.N(:))];
        
        hold(ax2, 'on')
        if max(ur)>0
            ylim(ax2,[ -max(ur) max(ur)])
        end
        %ylim(ax2,[ -max(ur) max(ur)])
        ylimits=ylim;
        plot(ax2,M.rgrid(id1)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        plot(ax2,M.rgrid(id2)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        
        ax3=subplot(2,2,3,'Parent',fig);
        uthet=M.fluidUTHET(:,zpos,fieldstep);
        plot(ax3,M.rgrid*1e3,uthet,'linewidth',1.5);
        xlim(ax3,[M.rgrid(1) M.rgrid(end)]*1e3)
        xlabel(ax3,'r [mm]')
        title(ax3,'Azimuthal velocity')
        ylabel(ax3,'v_\theta [m/s]');
        %c.Limits=[0 max(M.N(:))];
        
        hold(ax3, 'on')
        if max(uthet)>0
            ylim(ax3,[ -max(uthet) max(uthet)])
        end
        ylimits=ylim;
        plot(ax3,M.rgrid(id1)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        plot(ax3,M.rgrid(id2)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        
        uExb=-M.Er(:,zpos,fieldstep)./M.Bz(zpos,:)'.*(uthet~=0);
        plot(ax3,M.rgrid*1e3,uExb,'linewidth',1.5);
        
        ax4=subplot(2,2,4,'Parent',fig);
        uz=M.fluidUZ(:,zpos,fieldstep);
        plot(ax4,M.rgrid*1e3,uz,'linewidth',1.5);
        xlim(ax4,[M.rgrid(1) M.rgrid(end)]*1e3)
        xlabel(ax4,'r [mm]')
        title(ax4,'Axial velocity')
        ylabel(ax4,'v_z [m/s]');
        %c.Limits=[0 max(M.N(:))];
        
        hold(ax4, 'on')
        if max(uz)>0
            ylim(ax4,[ -max(uz) max(uz)])
        end
        ylimits=ylim;
        plot(ax4,M.rgrid(id1)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        plot(ax4,M.rgrid(id2)*[1 1]*1e3,ylimits,'k--','linewidth',1.5,'Displayname','Boundaries');
        
        linkaxes([ax1,ax2,ax3,ax4],'x');
    end

    function plotGridButtonPushed(btn,ax)
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        f=figure();
        PlotEspic2dgriddata(f,M,sld.Value,edr.Value);
        f.PaperOrientation='landscape';
        [~, name, ~] = fileparts(M.file);
        print(f,sprintf('%sGrid%d%d',name,sld.Value,edr.Value),'-dpdf','-fillpage')
    end

    function updatefigdata(control, event, Othercontrol, fig)
        
        if contains(event.Source.Tag,'time')
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
        elseif contains(event.Source.Tag,'axial')
            if strcmp(event.EventName,'ValueChanged')
                zpos=floor(control.Value);
                control.Value=zpos;
            elseif strcmp(event.EventName,'KeyPress')
                zpos=floor(control.Value);
                control.Value=zpos;
            else
                zpos=floor(event.Value);
            end
            Othercontrol.Value=zpos;
        end
        
        updatesubplotsdata(fieldstep, zpos, fig);
    end

    function updatesubplotsdata(fieldstep, zpos, fig)
        sgtitle(fig,sprintf('step=%d t=%0.5e s z=%0.3e mm',(fieldstep-1)*M.it1,M.t2d(fieldstep),M.zgrid(zpos)*1e3))
        [~,rcenterid]=max(M.geomweight(:,zpos,1));
        [~,id1]=min(abs(M.geomweight(1:rcenterid,zpos,1)));
        [~,id2]=min(abs(M.geomweight(rcenterid:end,zpos,1)));
        id2=id2+rcenterid;
        rlim1=M.rgrid(id1)*[1 1]*1e3;
        rlim2=M.rgrid(id2)*[1 1]*1e3;
        
        %% update density
        ax1=fig.Children(end);
        geomw=M.geomweight(:,zpos,1)>=0;
        dens=M.N(:,zpos,fieldstep).*geomw;
        Er=M.Er(:,zpos,fieldstep).*geomw;
        Ez=M.Ez(:,zpos,fieldstep).*geomw;
        yyaxis(ax1,'left')
        ax1.Children(end).YData=dens;
        ylimits=ylim(ax1);
        ax1.Children(end-1).XData=rlim1;
        ax1.Children(end-1).YData=ylimits;
        ax1.Children(end-2).XData=rlim2;
        ax1.Children(end-2).YData=ylimits;
        
        yyaxis(ax1,'right')
        ax1.Children(end).YData=Er;
        ax1.Children(end-1).YData=Ez;
        if max(abs([Er; Ez]))>0 
            ylim(ax1,max(abs([Er; Ez]))*[ -1 1])
        end
        
        %     view(ax1,2)
        %% update Radial velocity
        ax2=fig.Children(end-2);
        ur=M.fluidUR(:,zpos,fieldstep).*geomw;
        
        ax2.Children(end).YData=ur;
        if max(abs(ur))>0
            ylim(ax2,max(abs(ur))*[ -1 1])
        end
        ax2.Children(end-1).XData=rlim1;
        ax2.Children(end-2).XData=rlim2;
        ylimits=ylim(ax2);
        ax2.Children(end-1).YData=ylimits;
        ax2.Children(end-2).YData=ylimits;
        
        
        %% update Azimuthal velocity
        ax3=fig.Children(end-3);
        uthet=M.fluidUTHET(:,zpos,fieldstep).*geomw;
        uExb=-M.Er(:,zpos,fieldstep)./M.Bz(zpos,:)'.*(uthet~=0);
        
        ax3.Children(end-3).YData=uExb';
        ax3.Children(end).YData=uthet;
        
        if max(uthet)>0
            ylim(ax3,[ -max(uthet) max(uthet)])
        end
        
        ax3.Children(end-1).XData=rlim1;
        ax3.Children(end-2).XData=rlim2;
        ylimits=ylim(ax3);
        ax3.Children(end-1).YData=ylimits;
        ax3.Children(end-2).YData=ylimits;
        
       %% update Axial velocity
        ax4=fig.Children(end-4);
        uz=M.fluidUZ(:,zpos,fieldstep).*geomw;
        
        
        ax4.Children(end).YData=uz;
        
        if max(abs(uz))>0
            ylim(ax4,max(abs(uz))*[ -1 1])
        end
        ax4.Children(end-1).XData=rlim1;
        ax4.Children(end-2).XData=rlim2;
        ylimits=ylim(ax4);
        ax4.Children(end-1).YData=ylimits;
        ax4.Children(end-2).YData=ylimits;
        drawnow limitrate
        
    end

end


