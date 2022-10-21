function f=dispespicFluid(M,flux,log,parper)
fieldstep=length(M.t2d);
if nargin<2
    flux=false;
end
if nargin<3
    log=false;
end
if nargin<4
    parper=false;
end

f=uifigure('Name','Fluid data');

mf=uipanel(f,'Position',[5 50 f.Position(3)-10 f.Position(4)-55]);
mf.AutoResizeChildren='off';
m=uipanel(f,'Position',[5 5 f.Position(3)-10 40]);

sgtitle(mf,sprintf('t=%0.5e s',M.t2d(fieldstep)))

sld = uislider(m,'Position',[10 30 0.6*m.Position(3) 3]);
sld.Limits=[1 length(M.t2d)];
sld.Value=fieldstep;


edt = uieditfield(m,'numeric','Limits',[1 length(M.t2d)],'Value',1);
edt.Position=[sld.Position(1)+sld.Position(3)+25 5 40 20];
edt.RoundFractionalValues='on';


Printbt=uibutton(m,'Position',[edt.Position(1)+edt.Position(3)+10 5 40 20],'Text', 'Save');
%Playbt=uibutton(m,'Position',[Printbt.Position(1)+Printbt.Position(3)+10 5 40 20],'Text', 'Play/Pause');


sld.ValueChangingFcn={@updatefigdata,edt,mf};
edt.ValueChangedFcn={@updatefigdata,sld,mf};
Printbt.ButtonPushedFcn={@plotGridButtonPushed};

set(f,'KeyPressFcn',{ @onKeyDown,sld,edt,mf})

[R,Z]=meshgrid(M.rgrid,M.zgrid);
Rinv=1./R;
Rinv(:,1)=0;

PlotEspic2dfluiddata(mf,M,fieldstep);

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

function PlotEspic2dfluiddata(fig,M,fieldstep)
%PlotEspic2dgriddata Plot the 2d data of espic2d at time step fieldstep
sgtitle(fig,sprintf('t=%0.5e s',M.t2d(fieldstep)))


ax1=subplot(2,2,1,'Parent',fig);
N=M.N(:,:,fieldstep);
N(M.geomweight(:,:,1)<0)=NaN;
contourf(ax1,M.zgrid,M.rgrid,N,50,'linecolor','none')
xlim(ax1,[M.zgrid(1) M.zgrid(end)])
ylim(ax1,[M.rgrid(1) M.rgrid(end)])
xlabel(ax1,'z [m]')
ylabel(ax1,'r [m]')
title(ax1,'Density')
c = colorbar(ax1);
c.Label.String= 'n[m^{-3}]';
%c.Limits=[0 max(M.N(:))];
%caxis(ax1,[0 max(M.N(:))]);
view(ax1,2)

%set(ax1,'colorscale','log')

UR0=M.fluidUR(:,:,fieldstep);
UZ0=M.fluidUZ(:,:,fieldstep);

if parper
    UR=UR0.*M.sinthet-UZ0.*M.costhet;    
    UZ=UR0.*M.costhet+UZ0.*M.sinthet;
else
    UR=UR0;
    UZ=UZ0;
end

UR(M.geomweight(:,:,1)<0)=NaN;
ax2=subplot(2,2,2,'Parent',fig);
if flux
    UR=UR.*M.N(:,:,fieldstep);
end
if log
    UR=abs(UR);
end
contourf(ax2,M.zgrid,M.rgrid,UR,50,'linecolor','none');
%plot(ax2,M.zgrid,data.pot(:,5))
xlabel(ax2,'z [m]')
ylabel(ax2,'r [m]')
colormap(ax2,'jet')
c = colorbar(ax2);
%c.Limits=[min(M.fluidUR(:,:,:)) max(M.fluidUR(:,:,:))];
titl='';
labl='';
if parper
    labl='U_{per}';
    titl= 'Perpendicular';
else
    labl='U_{r}';
    titl= 'Radial';
end
if flux
    labl= [labl ' [1/sm^2]'];
    titl=[titl  ' Flux'];
    caxis(ax2,5e21*[-1 1])
else
    labl= [labl ' [m/s]'];
    titl=[titl  ' velocity'];
    %caxis(ax2,[-4e5 2e4])
end
if log
    labl= ['|' labl '|'];
    set(ax2,'ColorScale','log') 
end
title(ax2,titl)
c.Label.String=labl;
grid(ax2, 'on')
%caxis(ax2,[-4e5 1e6])

view(ax2,2)

ax3=subplot(2,2,3,'Parent',fig);
omegathet=M.fluidUTHET(:,:,fieldstep);
if flux
    omegathet=omegathet.*M.N(:,:,fieldstep);
else
    %omegathet=abs(omegathet.*Rinv');
end
omegathet(M.geomweight(:,:,1)<0)=NaN;
contourf(ax3,M.zgrid,M.rgrid,omegathet,50,'linecolor','none')

%mean(M.fluidUTHET(:,:,end).*Rinv')
xlabel(ax3,'z [m]')
ylabel(ax3,'r [m]')
colormap(ax3,'jet')
c = colorbar(ax3);
if flux
    c.Label.String= 'U_\theta [1/sm^2]';
    title(ax3,'Azimuthal Flux')
else
    c.Label.String= '\omega_\theta [m/s]';
    title(ax3,'Azimuthal rotation velocity')
    
end

if log
    set(ax3,'colorscale','log')
    end

grid(ax3, 'on')
view(ax3,2)

ax4=subplot(2,2,4,'Parent',fig);
UZ(M.geomweight(:,:,1)<0)=NaN;
if flux
    UZ=UZ.*M.N(:,:,fieldstep);
end
if log
    UZ=abs(UZ);
end
contourf(ax4,M.zgrid,M.rgrid,UZ,30,'linecolor','none')
xlabel(ax4,'z [m]')
ylabel(ax4,'r [m]')
colormap(ax4,'jet')
c = colorbar(ax4);
%c.Limits=[min(M.fluidUZ(:)) max(M.fluidUZ(:))];
titl='';
labl='';
if parper
    labl='U_{par}';
    titl= 'Parallel';
else
    labl='U_{z}';
    titl= 'Axial';
end
if flux
    labl= [labl ' [1/sm^2]'];
    titl=[titl  ' Flux'];
else
    labl= [labl ' [m/s]'];
    titl=[titl  ' velocity'];
end
if log
    labl= ['|' labl '|'];
    set(ax2,'ColorScale','log') 
end
title(ax4,titl)
c.Label.String=labl;
grid(ax4, 'on')
%caxis(ax4,[min(M.fluidUZ(:)) max(M.fluidUZ(:))])
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
    elseif strcmp(event.EventName,'KeyPress')
       fieldstep=floor(control.Value);
       control.Value=fieldstep;
    else
      fieldstep=floor(event.Value);
    end
    Othercontrol.Value=fieldstep;
    
    sgtitle(fig,sprintf('t=%0.5e s',double((fieldstep-1)*M.it1)*M.dt))
    
    %% update Position histogram
    ax1=fig.Children(end);
    N=M.N(:,:,fieldstep);
    N(M.geomweight(:,:,1)<0)=NaN;
    %ax1.Children(1).CData=N;
    ax1.Children(1).ZData=N;
    
    
    UR0=M.fluidUR(:,:,fieldstep);
    UZ0=M.fluidUZ(:,:,fieldstep);
    if parper
    UR=UR0.*M.sinthet-UZ0.*M.costhet;    
    UZ=UR0.*M.costhet+UZ0.*M.sinthet;
else
    UR=UR0;
    UZ=UZ0;
    end
    UR(M.geomweight(:,:,1)<0)=NaN;
    UZ(M.geomweight(:,:,1)<0)=NaN;
    view(ax1,2)
    %% update Radial velocity

if flux
    UR=UR.*N;
end
if log
    UR=abs(UR);
end
    %fig.Children(end-2).Children(1).CData=UR;
    fig.Children(end-2).Children(1).ZData=UR;
    %% update Azimuthal velocity
omegathet=M.fluidUTHET(:,:,fieldstep);
omegathet(M.geomweight(:,:,1)<0)=NaN;
if flux
    omegathet=omegathet.*M.N(:,:,fieldstep);
else
    %omegathet=abs(omegathet.*Rinv');
end
    %fig.Children(end-4).Children(1).CData=omegathet;
    fig.Children(end-4).Children(1).ZData=omegathet;
    %% update Axial velocity
if flux
    UZ=UZ.*N;
end
if log
    UZ=abs(UZ);
end
    %fig.Children(end-6).Children(1).CData=UZ;
    fig.Children(end-6).Children(1).ZData=UZ;
    %drawnow limitrate
    
    

end

end


