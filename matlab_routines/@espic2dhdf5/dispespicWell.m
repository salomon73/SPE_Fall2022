function dispespicWell(M,fracn,logdensity,showgrid,fixed)
%dispespicFields Allows to display the time evolution of the density, and Penning potential well
%   M is of class espic2dhdf5 and contains the simulation results
fieldstep=1;
if nargin <3
    logdensity=false;
end
if nargin <4
    showgrid=false;
end
if nargin <5
    fixed=false;
end
if nargin<2
fracn=0.1;
end

fixed=fi(fixed);

f=uifigure('Name',sprintf('Well data %s',M.name));

mf=uipanel(f,'Position',[5 50 f.Position(3)-10 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);

sgtitle(mf,sprintf('step=%d t=%0.5e s',fieldstep*M.it1,M.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.4*m.Position(3) 3]);
sld.Value=fieldstep;
sld.Limits=[1 size(M.t2d,1)];

edt = uieditfield(m,'numeric','Limits',[1 size(M.t2d,1)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';

MaxN=0;
Minwell=0;
Maxwell=0;
plotaxes=gobjects(2);

Printbt=uibutton(m,'Position',[edt.Position(1)+edt.Position(3)+10 5 40 20],'Text', 'Save');
Play=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play');
Pause=uibutton(m,'Position',[Play.Position(1)+Play.Position(3)+10 5 45 20],'Text', 'Pause');
rcoordbt=uiswitch(m);
rcoordbt.Items = {'rz','rpsi'};
rcoordbt.Position(1:2)=[Pause.Position(1)+Pause.Position(3)+25 5];

stop=false;

sld.ValueChangingFcn={@updatefigdata,edt,mf};
edt.ValueChangedFcn={@updatefigdata,sld,mf};
Printbt.ButtonPushedFcn={@plotGridButtonPushed};
Play.ButtonPushedFcn={@plotPlayButtonPushed};
rcoordbt.ValueChangedFcn={@change_radial_coord};
rcoordbt.Value='rz';

Pause.ButtonPushedFcn={@PauseButtonPushed};
PlotEspic2dgriddata(mf,M,fieldstep);

    function plotPlayButtonPushed(btn,ax)
        stop=false;
        for i=edt.Value:5:sld.Limits(2)
            edt.Value=i;
            sld.Value=i;
            updatesubplotsdata(i,mf);
            pause(0.01)
            if stop
                stop=false;
                break;
            end
            
        end
    end
    function PauseButtonPushed(btn,ax)
        stop = true;
    end
    function change_radial_coord(btn,ax)
        if strcmp(rcoordbt.Value,'rz')
            ylabel(plotaxes(2),'r [m]')
            ylim(plotaxes(2),[M.rgrid(1) M.rgrid(end)])
            ylabel(plotaxes(1),'r [m]')
            ylim(plotaxes(1),[M.rgrid(1) M.rgrid(end)])
            linkaxes(plotaxes)
        else
            ylabel(plotaxes(2),'rA_\theta [Tm^2]')
            ylim(plotaxes(2),[M.rAthet(1,1) M.rAthet(end,1)])
            ylabel(plotaxes(1),'rA_\theta [Tm^2]')
            ylim(plotaxes(1),[M.rAthet(1,1) M.rAthet(end,1)])
            %linkaxes(plotaxes,'off');
            %ylim(plotaxes(1),[M.rgrid(1) M.rgrid(end)])
        end
        updatesubplotsdata(sld.Value, mf);
    end

    function PlotEspic2dgriddata(fig,M,fieldstep)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*M.it1,M.t2d(fieldstep)))
        
        N=M.N(:,:,fieldstep);
        plotaxes(1)=subplot(2,1,1,'Parent',fig);
        sf=surface(plotaxes(1),M.zgrid,M.rgrid,N,'edgecolor','none');
        if(showgrid)
            set(sf,'edgecolor','black')
        end
        hold(plotaxes(1), 'on')
        contour(plotaxes(1),M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r--')
        xlim(plotaxes(1),[M.zgrid(1) M.zgrid(end)])
        ylim(plotaxes(1),[M.rgrid(1) M.rgrid(end)])
        xlabel(plotaxes(1),'z [m]')
        ylabel(plotaxes(1),'r [m]')
        title(plotaxes(1),'Density')
        c = colorbar(plotaxes(1));
        c.Label.String= 'n[m^{-3}]';
        %c.Limits=[0 max(M.N(:))];
        if(isboolean(fixed) && fixed)
            climits=caxis(plotaxes(1));
            MaxN=climits(2);
        elseif(~isboolean(fixed))
            MaxN=fixed.data;
            caxis(plotaxes(1),[-Inf MaxN]);
        end
        view(plotaxes(1),2)
        
        if logdensity
            set(plotaxes(1),'zscale','log')
            set(plotaxes(1),'colorscale','log')
        end
        
        
        plotaxes(2)=subplot(2,1,2,'Parent',fig);
        model=M.potentialwellmodel(fieldstep);
        z=model.z;
        r=model.r;
        pot=model.pot;
        rathet=model.rathet;
        dispz=M.zgrid;
        if rcoordbt.Value
            [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rgrid);
            pot=griddata(z,r,pot,Zmesh,Rmesh);
            dispr=M.rgrid;
        else
            [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,1));
            pot=griddata(z,rathet,pot,Zmesh,Rmesh);
            N=griddata(Zmesh,M.rAthet,N,Zmesh,Rmesh);
            dispr=M.rAthet(:,1);
        end
        
        sf=contourf(plotaxes(2),dispz,dispr,pot(1:end,1:end),'edgecolor','none');
        hold(plotaxes(2), 'on')
        contour(plotaxes(2),dispz,dispr,N-fracn*max(N(:)),[0 0],'k--')
        contour(plotaxes(2),dispz,dispr,M.geomweight(:,:,1),[0 0],'r--')
        if(showgrid)
            set(sf,'edgecolor','black')
        end
        xlim(plotaxes(2),[M.zgrid(1) M.zgrid(end)])
        xlabel(plotaxes(2),'z [m]')
        if rcoordbt.Value
            ylabel(plotaxes(2),'r [m]')
            ylim(plotaxes(2),[M.rgrid(1) M.rgrid(end)])
        else
            ylabel(plotaxes(2),'rA_\theta [Tm^2]')
            ylim(plotaxes(2),[M.rAthet(1,1) M.rAthet(end,1)])
        end
        title(plotaxes(2),'Well')
        c = colorbar(plotaxes(2));
        c.Label.String= 'depth [eV]';
        %c.Limits=[0 max(M.N(:))];
        if(isboolean(fixed) && fixed)
            climits=caxis(plotaxes(2));
            Minwell=climits(1);
            Maxwell=climits(2);
        elseif(~isboolean(fixed))
            climits=caxis(plotaxes(2));
            Minwell=climits(1);
            Maxwell=climits(2);
        end
        view(plotaxes(2),2)
        colormap(plotaxes(2),'jet');
        Blines=M.rAthet;
        levels=linspace(min(Blines(M.geomweight(:,:,1)>0)),max(Blines(M.geomweight(:,:,1)>0)),20);
        Blines(M.geomweight(:,:,1)<0)=NaN;
        contour(plotaxes(2),M.zgrid,M.rgrid,Blines,real(levels),'m-.','linewidth',1.5,'Displayname','Magnetic field lines');
        
        linkaxes(plotaxes);
        axis(plotaxes(1), 'equal')
        axis(plotaxes(2), 'equal')
    end

    function plotGridButtonPushed(btn,ax)
        f=figure();
        PlotEspic2dgriddata(f,M,edt.Value);
        f.PaperOrientation='landscape';
        [~, name, ~] = fileparts(M.file);
        f.PaperSize=[18,12];
        savefig(f,sprintf('%sGridFields%d',name,floor(edt.Value)))
        print(f,sprintf('%sGridWell%d',name,edt.Value),'-dpdf','-fillpage')
    end

    function updatefigdata(control, event, Othercontrol, fig)
        
        
        if strcmp(event.EventName,'ValueChanged')
            fieldstep=floor(control.Value);
            control.Value=fieldstep;
        else
            fieldstep=floor(event.Value);
        end
        Othercontrol.Value=fieldstep;
        updatesubplotsdata(fieldstep, fig);
    end

    function updatesubplotsdata(fieldstep, fig)
        
        sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*M.it1,double((fieldstep-1)*M.it1)*M.dt))
        
        %% update Position histogram
        dens=M.N(:,:,fieldstep);
        
        model=M.potentialwellmodel(fieldstep);
        z=model.z;
        r=model.r;
        pot=model.pot;
        rathet=model.rathet;
        geomweight=M.geomweight(:,:,1);
        if strcmp(rcoordbt.Value,'rz')
            plotaxes(1).Children(end).YData=M.rgrid(1:end);
            gweight=geomweight;
            axis(plotaxes(1),'equal')
        else
            [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,1));
            dens=griddata(Zmesh(:),M.rAthet(:),dens(:),Zmesh,Rmesh);
            gweight=griddata(Zmesh(:),M.rAthet(:),geomweight(:),Zmesh,Rmesh);
            plotaxes(1).Children(end).YData=M.rAthet(:,1);
            axis(plotaxes(1),'normal')
        end
        
        dens(gweight<0)=NaN;
        plotaxes(1).Children(end).ZData=dens;
        plotaxes(1).Children(end).CData=dens;
        if(isboolean(fixed) && fixed)
            climits=caxis(plotaxes(1));
            MaxN=climits(2);
            nmax=max(dens(:));
            if(nmax>MaxN)
                MaxN=nmax;
            end
            caxis(plotaxes(1),[0 MaxN]);
        elseif(~isboolean(fixed))
            MaxN=fixed.data;
            caxis(plotaxes(1),[-Inf MaxN]);
        end
        
        if strcmp(rcoordbt.Value,'rz')
            [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rgrid);
            well=griddata(z,r,pot,Zmesh,Rmesh);
            
            plotaxes(2).Children(end).YData=M.rgrid(1:end);
            plotaxes(2).Children(end-1).YData=M.rgrid(1:end);
            axis(plotaxes(2),'equal')
        else
            [Zmesh,Rmesh]=meshgrid(M.zgrid,M.rAthet(:,1));
            well=griddata(z,rathet,pot,Zmesh,Rmesh);
            
            dens=griddata(Zmesh,M.rAthet,dens,Zmesh,Rmesh);
            
            plotaxes(2).Children(end).YData=M.rAthet(:,1);
            plotaxes(2).Children(end-1).YData=M.rAthet(:,1);
            axis(plotaxes(2),'normal')
        end
        well(gweight<0)=NaN;
        
        plotaxes(2).Children(end-1).ZData=dens-fracn*max(dens(:));
        
        plotaxes(2).Children(end).ZData=well;
        %plotaxes(2).Children(end).CData=well;
        if(isboolean(fixed) && fixed || ~isboolean(fixed))
            climits=caxis(plotaxes(2));
            Maxwell=max([well(:);climits(2)]);
            Minwell=min([well(:);climits(1)]);
            caxis(plotaxes(2),[Minwell Maxwell]);
        end
        drawnow limitrate
        
    end

end


