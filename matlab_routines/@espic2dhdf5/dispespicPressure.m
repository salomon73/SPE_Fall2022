function dispespicPressure(M,logdensity,showgrid,fixed,temperature)
%dispespicPressure Allows to display the time evolution of the pressure tensor
%   M is of class espic2dhdf5 and contains the simulation results
fieldstep=1;
if nargin <2
    logdensity=false;
end
if nargin <3
    showgrid=false;
end
if nargin <4
    fixed=false;
end
if nargin <5
    temperature=false;
end
fixed=fi(fixed);

f=uifigure('Name',sprintf('Pressure data %s',M.name));

mf=uipanel(f,'Position',[5 50 f.Position(3)-10 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);

sgtitle(mf,sprintf('step=%d t=%0.5e s',fieldstep*M.it1,M.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.6*m.Position(3) 3]);
sld.Value=fieldstep;
sld.Limits=[1 size(M.t2d,1)];

edt = uieditfield(m,'numeric','Limits',[1 size(M.t2d,1)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';

MaxP=0;
MinP=0;
Paxes=gobjects(6);

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
PlotEspic2dgriddata(mf,M,fieldstep);

    function plotPlayButtonPushed(btn,ax)
        stop=false;
        for i=edt.Value:sld.Limits(2)
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

    function PlotEspic2dgriddata(fig,M,fieldstep)
        %PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
        sgtitle(fig,sprintf('step=%d t=%0.5e s',(fieldstep-1)*M.it1,M.t2d(fieldstep)))
        if temperature
            Clabels={'T_{rr} [eV]', 'T_{r\theta} [eV]', 'T_{rz} [eV]', 'T_{\theta\theta} [eV]', 'T_{\theta z} [eV]', 'T_{zz} [eV]'};
        else
            Clabels={'P_{rr} [Pa]', 'P_{r\theta} [Pa]', 'P_{rz} [Pa]', 'P_{\theta\theta} [Pa]', 'P_{\theta z} [Pa]', 'P_{zz} [Pa]'};
        end
        for i=1:6
            Paxes(i)=subplot(2,3,i,'Parent',fig);
            if temperature
                N=M.N(:,:,fieldstep);
                invN=1./N;
                invN(isinf(invN))=0;
                p=squeeze(M.Presstens(i,:,:,fieldstep)).*invN/M.qe;
            else
                p=squeeze(M.Presstens(i,:,:,fieldstep));
            end
            if logdensity
                p=abs(p);
            end
            contourf(Paxes(i),M.zgrid,M.rgrid,p,'edgecolor','none');
            hold(Paxes(i),'on')
            if(showgrid)
                set(sf,'edgecolor','black')
            end
            contour(Paxes(i),M.zgrid,M.rgrid,M.geomweight(:,:,1),[1e-8 1e-8],'r-','linewidth',2.5,'Displayname','Boundaries');
            xlim(Paxes(i),[M.zgrid(1) M.zgrid(end)])
            ylim(Paxes(i),[M.rgrid(1) M.rgrid(end)])
            xlabel(Paxes(i),'z [m]')
            ylabel(Paxes(i),'r [m]')
            title(Paxes(i),Clabels{i})
            c = colorbar(Paxes(i));
            c.Label.String= Clabels{i};
            if(isboolean(fixed) && fixed)
                climits=caxis(Paxes(i));
                MaxP=max(MaxP,climits(2));
                MinP=min(MinP,climits(1));
            elseif(~isboolean(fixed))
                MaxP=fixed.data;
                MinP=-Inf;
                caxis(ax1,[MinP MaxP]);
            end
            
            if logdensity
                %set(Paxes(i),'zscale','log')
                set(Paxes(i),'colorscale','log')
            end
            colormap(Paxes(i),cool)
        end
        if(isboolean(fixed) && fixed)
            for i=1:6
                caxis(Paxes(i),[MinP MaxP]);
            end
        end
        linkaxes(Paxes);
        
    end

    function plotGridButtonPushed(btn,ax)
        %UNTITLED2 Summary of this function goes here
        %   Detailed explanation goes here
        f=figure();
        PlotEspic2dgriddata(f,M,sld.Value);
        f.PaperOrientation='landscape';
        [~, name, ~] = fileparts(M.file);
        print(f,sprintf('%sGrid%d',name,sld.Value),'-dpdf','-fillpage')
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
        
        %% update Pressure subplots
        for i=1:6
            if temperature
                N=M.N(:,:,fieldstep);
                invN=1./N;
                invN(isinf(invN))=0;
                P=squeeze(M.Presstens(i,:,:,fieldstep)).*invN/M.qe;
            else
                P=squeeze(M.Presstens(i,:,:,fieldstep));
            end
            Paxes(i).Children(end).ZData=P;
            %Paxes(i).Children.CData=P;
            if(isboolean(fixed) && fixed)
                climits=caxis(Paxes(i));
                MaxP=max([P(:);climits(2)]);
                MinP=min([P(:);climits(1)]);
                caxis(Paxes(i),[0 MaxP]);
            elseif(~isboolean(fixed))
                MaxP=fixed.data;
                caxis(Paxes(i),[-Inf MaxP]);
            end
        end
        if(isboolean(fixed) && fixed)
            for i=1:6
                caxis(Paxes(i),[MinP MaxP]);
            end
        end
        drawnow limitrate
        
    end

end


