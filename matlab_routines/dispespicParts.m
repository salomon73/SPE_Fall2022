function dispespicParts(M,fieldstep,parper)
%dispespicParts Show on an interactive plot the phase space of the given specie in M
%   Detailed explanation goes here

if nargin<2
    fieldstep=1;
end
if nargin<3
    parper=false;
end
f2=uifigure('Name','Particles data');

mf2=uipanel(f2,'Position',[5 50 f2.Position(3)-10 f2.Position(4)-55]);
mf2.AutoResizeChildren='off';
m2=uipanel(f2,'Position',[5 5 f2.Position(3)-10 40]);

sgtitle(mf2,sprintf('step=%d t=%0.5e s',fieldstep*M.it2,M.tpart(fieldstep)))

sld2 = uislider(m2,'Position',[10 30 0.6*m2.Position(3) 3]);
sld2.Value=fieldstep;
sld2.Limits=[1 length(M.tpart)];

edt2 = uieditfield(m2,'numeric','Limits',[1 length(M.tpart)],'Value',1);
edt2.Position=[sld2.Position(1)+sld2.Position(3)+25 5 40 20];
edt2.RoundFractionalValues='on';

Printbt2=uibutton(m2,'Position',[edt2.Position(1)+edt2.Position(3)+10 5 40 20],'Text', 'Save');

sld2.ValueChangingFcn={@updatefigpartdata,edt2,mf2};
edt2.ValueChangedFcn={@updatefigpartdata,sld2,mf2};

Printbt2.ButtonPushedFcn={@plotPartButtonPushed};

set(f2,'KeyPressFcn',{ @onKeyDown,sld2,edt2,mf2})


PlotEspic2dpartdata(mf2,M,fieldstep);

    function PlotEspic2dpartdata(fig,M,fieldstep)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        %sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*M.it2,M.tpart(fieldstep)))
        sgtitle(fig,sprintf('t=%0.5e s nbparts=%d',M.tpart(fieldstep),M.nbparts(fieldstep)))
        
        vmax=4e7;
        
        Z=M.Z(1:min(M.nbparts(fieldstep),M.Z.nparts),fieldstep);
        R=M.R(1:min(M.nbparts(fieldstep),M.R.nparts),fieldstep);
        
        if parper
            vpar=M.Vpar(1:min(M.nbparts(fieldstep),M.Z.nparts),fieldstep);
            vper=M.Vperp(1:min(M.nbparts(fieldstep),M.Z.nparts),fieldstep);
        else
            vpar=M.VZ(1:min(M.nbparts(fieldstep),M.Z.nparts),fieldstep);
            vper=M.VR(1:min(M.nbparts(fieldstep),M.Z.nparts),fieldstep);
        end
        
        ax1=subplot(3,2,1,'Parent',fig);
        if isa(M,'h5parts')
            contour(ax1,M.zgrid,M.rgrid,M.parent.geomweight(:,:,1),[0 0],'r--')
            rindex=M.rindex;
            zindex=M.zindex;
        else
            contour(ax1,M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0],'r--')
            rindex=[1 length(M.rgrid)];
            zindex=[1 length(M.zgrid)];
        end
        hold(ax1, 'on')
        rectangle(ax1,'Position',[M.zgrid(zindex(1)) M.rgrid(rindex(1))...
            M.zgrid(zindex(end))-M.zgrid(zindex(1)) M.rgrid(rindex(end))-M.rgrid(rindex(1))],'linestyle','--')
        plot(ax1,Z,R,'.');
        
        if isa(M,'h5parts')
        xlim(ax1,[M.zgrid(zindex(1)) M.zgrid(zindex(end))])
        ylim(ax1,[M.rgrid(rindex(1)) M.rgrid(rindex(end))])
        end
        xlabel(ax1,'Z [m]')
        ylabel(ax1,'R [m]')
        title(ax1,'Position')
        grid(ax1, 'on')
        
        ax2=subplot(3,2,2,'Parent',fig);
        plot(ax2,vpar,vper,'.');
        if parper
            xlabel(ax2,'V_{par} [m/s]')
            ylabel(ax2,'V_{perp} [m/s]')
        else
            xlabel(ax2,'V_Z[m/s]')
            ylabel(ax2,'V_R [m/s]')
        end
        ylim(ax2,vmax*[-1 1])
        xlim(ax2,vmax*[-1 1])
        axis(ax2, 'equal')
        grid(ax2, 'on')
        
        
        ax3=subplot(3,2,3,'Parent',fig);
        plot(ax3,Z,vpar,'.');
        xlabel(ax3,'Z [m]')
        if parper
            ylabel(ax3,'V_{par} [m/s]')
        else
            ylabel(ax3,'V_Z[m/s]')
        end
        grid(ax3, 'on')
        xlim(ax3,[M.zgrid(zindex(1)) M.zgrid(zindex(end))])
        ylim(ax3,vmax*[-1 1])
        
        ax4=subplot(3,2,4,'Parent',fig);
        plot(ax4,R,vpar,'.')
        if parper
            ylabel(ax4,'V_{par} [m/s]')
        else
            ylabel(ax4,'V_Z[m/s]')
        end
        xlabel(ax4,'R [m]')
        xlim(ax4,[M.rgrid(rindex(1)) M.rgrid(rindex(end))])
        ylim(ax4,vmax*[-1 1])
        grid(ax4, 'on')
        
        ax5=subplot(3,2,5,'Parent',fig);
        plot(ax5,Z,vper,'.');
        xlim(ax5,[M.zgrid(zindex(1)) M.zgrid(zindex(end))])
        xlabel(ax5,'Z [m]')
        if parper
            ylabel(ax5,'V_{perp} [m/s]')
        else
            ylabel(ax5,'V_R [m/s]')
        end
        grid(ax5, 'on')
        ylim(ax5,vmax*[-1 1])
        
        ax6=subplot(3,2,6,'Parent',fig);
        plot(ax6,R,vper,'.');
        xlim(ax6,[M.rgrid(rindex(1)) M.rgrid(rindex(end))])
        if parper
            ylabel(ax6,'V_{perp} [m/s]')
        else
            ylabel(ax6,'V_R [m/s]')
        end
        xlabel(ax6,'R [m]')
        ylim(ax6,vmax*[-1 1])
        grid(ax6, 'on')
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
            updatefigpartdata(slider, event, editfield ,fig)
        end
    end


    function updatefigpartdata(control, event, Othercontrol, fig)
        
        
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
        
        sgtitle(fig,sprintf('t=%0.5e s nbparts=%03d',M.tpart(fieldstep),M.nbparts(fieldstep)))
        
        %         if isa(M,'espic2dhdf5')
        %             [~,t0dstep]=min(abs(M.tpart(fieldstep)-M.t0d));
        %         else
        t0dstep=fieldstep;
        %         end
        
        Z=M.Z(1:min(M.nbparts(t0dstep),M.Z.nparts),fieldstep);
        R=M.R(1:min(M.nbparts(t0dstep),M.R.nparts),fieldstep);
        
        if parper
            vpar=M.Vpar(1:min(M.nbparts(t0dstep),M.Z.nparts),fieldstep);
            vper=M.Vperp(1:min(M.nbparts(t0dstep),M.Z.nparts),fieldstep);
        else
            vpar=M.VZ(1:min(M.nbparts(t0dstep),M.Z.nparts),fieldstep);
            vper=M.VR(1:min(M.nbparts(t0dstep),M.Z.nparts),fieldstep);
        end
        %% update Position plot
        ax1=fig.Children(end);
        ax1.Children(1).XData=Z;
        ax1.Children(1).YData=R;
        
        
        view(ax1,2)
        %% update VPAR VPERP
        fig.Children(end-1).Children(1).XData=vpar;
        fig.Children(end-1).Children(1).YData=vper;
        axis(fig.Children(end-1),'equal')
        
        %% update Z VZ
        fig.Children(end-2).Children(1).XData=Z;
        fig.Children(end-2).Children(1).YData=vpar;
        %% update VZ VR
        fig.Children(end-3).Children(1).XData=R;
        fig.Children(end-3).Children(1).YData=vpar;
        
        %% update R VR
        fig.Children(end-4).Children(1).XData=Z;
        fig.Children(end-4).Children(1).YData=vper;
        
        %% update Z VR
        fig.Children(end-5).Children(1).XData=R;
        fig.Children(end-5).Children(1).YData=vper;
        %drawnow limitrate
    end

    function plotPartButtonPushed(btn,ax)
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        f=figure();
        PlotEspic2dpartdata(f,M,sld.Value);
        f.PaperOrientation='portrait';
        [~, name, ~] = fileparts(M.file);
        print(f,sprintf('%sParts%d',name,sld.Value),'-dpdf','-fillpage')
    end

end

