function [Forces, Continuity]=FBspline_article(M,Forces, Continuity, Density,norm,log,parper,contourscale, zlims, rlims, clims)
%ForceBalance Show the radial force balance
%   Plot the three radial forces for the given time-step it or averaged
%   over the range of time steps defined in it

it=Forces.it;
switch nargin
    case 4
        parper=false;
        contourscale=0.2;
        norm=false;
        log=false;
        zlims=[-inf inf];
    case 5
        parper=false;
        log=false;
        contourscale=0.2;
        zlims=[-inf inf];
    case 6
        parper=false;
        contourscale=0.2;
        zlims=[-inf inf];
    case 7
        contourscale=0.2;
        zlims=[-inf inf];
    case 8
        zlims=[-inf inf];
    case 9
        rlims=[];
    case 10
end
if nargin<10
    rlims=[];
end

if nargin<11
    climsr=[-230000 190000];
    climsz=[-15000 13000];
    climsthet=[-1100 1800];
else
    climsr=clims(1,:);
    climsthet=clims(2,:);
    climsz=clims(3,:);
end

folder='Article_Forces';

[status, msg, msgID] = mkdir(folder);

contourcolor=[0 0 0];%[255 20 147]/255;
N=Density.N;
n=mean(N,3);

maxdens=max(N(:));

densitycontour=contourc(M.zgrid,M.rgrid,mean(N,3),[contourscale contourscale]*maxdens);

geometriccontour=contourc(M.zgrid,M.rgrid,M.geomweight(:,:,1),[0 0]);

Fmaxr=max([max(Forces.Eforcer(:)),max(Forces.Bforcer(:)),max(Forces.inertforcer(:)),max(Forces.Pressforcer(:))]);
Fmaxz=max([max(Forces.Eforcez(:)),max(Forces.Bforcez(:)),max(Forces.Pressforcez(:)),max(Forces.inertforcez(:))]);
FEr=mean(Forces.Eforcer,3);
FBr=mean(Forces.Bforcer,3);
FPr=mean(Forces.Pressforcer,3);
FIr=mean(Forces.inertforcer,3);
Fdr=mean(Forces.Dragforcer,3);

FBthet=mean(Forces.Bforcethet,3);
%FPthet=zeros(size(FBthet));%mean(Forces.Pressforcethet,3);
FPthet=mean(Forces.Pressforcethet,3);
FIthet=mean(Forces.inertforcethet,3);
Fdthet=mean(Forces.Dragforcethet,3);

FEz=mean(Forces.Eforcez,3);
FBz=mean(Forces.Bforcez,3);
FPz=mean(Forces.Pressforcez,3);
FIz=mean(Forces.inertforcez,3);
Fdz=mean(Forces.Dragforcez,3);

durdt=mean(Forces.durdt,3);
duthetdt=mean(Forces.duthetdt,3);
duzdt=mean(Forces.duzdt,3);

totforcer=(FEr+FBr+FPr+FIr+Fdr-durdt);
totforcethet=(FBthet+FPthet+FIthet+Fdthet-duthetdt);
totforcez=(FEz+FBz+FPz+FIz+Fdz-duzdt);
if log
    Er=Forces.Eforcer(Forces.Eforcer~=0);
    Br=Forces.Bforcer(Forces.Bforcer~=0);
    inertr=Forces.inertforcer(Forces.inertforcer~=0);
    pressr=Forces.Pressforcer(Forces.Pressforcer~=0);
    Fminr=min([min(abs(Er(:))),min(abs(Br(:))),min(abs(inertr(:))),min(abs(pressr(:)))]);
else
    Fminr=min([min(Forces.Eforcer(:)),min(Forces.Bforcer(:)),min(Forces.inertforcer(:)),min(Forces.Pressforcer(:))]);
end
rgridmax=sum(M.nnr(1:2));

if parper
    costhet=(M.Br./M.B)';
    sinthet=(M.Bz./M.B)';
    FErd=FEr.*sinthet-FEz.*costhet;
    FBrd=FBr.*sinthet-FBz.*costhet;
    FPrd=FPr.*sinthet-FPz.*costhet;
    FIrd=FIr.*sinthet-FIz.*costhet;
    Fdrd=Fdr.*sinthet-Fdz.*costhet;
    durd=durdt.*sinthet-duzdt.*costhet;
    
    FBthetd=FBthet;
    FPthetd=FPthet;
    FIthetd=FIthet;
    Fdthetd=Fdthet;
    duthetd=duthetdt;
    
    FEzd=FEr.*costhet+FEz.*sinthet;
    FBzd=FBr.*costhet+FBz.*sinthet;
    FPzd=FPr.*costhet+FPz.*sinthet;
    FIzd=FIr.*costhet+FIz.*sinthet;
    Fdzd=Fdr.*costhet+Fdz.*sinthet;
    duzd=durdt.*costhet+duzdt.*sinthet;
    
    totforcerd=totforcer.*sinthet-totforcez.*costhet;
    totforcethetd=totforcethet;
    totforcezd=totforcer.*costhet+totforcez.*sinthet;
else
    FErd=FEr;
    FBrd=FBr;
    FPrd=FPr;
    FIrd=FIr;
    Fdrd=Fdr;
    durd=durdt;
    
    FBthetd=FBthet;
    FPthetd=FPthet;
    FIthetd=FIthet;
    Fdthetd=Fdthet;
    duthetd=duthetdt;
    
    FEzd=FEz;
    FBzd=FBz;
    FPzd=FPz;
    FIzd=FIz;
    Fdzd=Fdz;
    duzd=duzdt;
    
    totforcerd=totforcer;
    totforcethetd=totforcethet;
    totforcezd=totforcez;
end

papsize=[12 17];

%% Radial forces
f=figure('Name','radial');


ax0=subplot(3,2,1);
plotvalue(ax0,durd,'mn\partial_tu_{r}', 'F_{a,r} [N/m^3]',true,climsr,zlims,rlims);

ax1=subplot(3,2,2);
plotvalue(ax1,FErd,'qnE_r', 'F_{E,r} [N/m^3]',true,climsr,zlims,rlims);

ax2=subplot(3,2,3);
plotvalue(ax2,FBrd,'qn(u_\thetaB_z)','F_{B,r} [N/m^3]',true,climsr,zlims,rlims);

ax3=subplot(3,2,4);
plotvalue(ax3,FIrd,'-mn(\bf{u\bullet\nabla})u_r','F_{i,r} [N/m^3]',true,climsr,zlims,rlims);

ax4=subplot(3,2,5);
plotvalue(ax4,FPrd,'-(\nablaP)_r','F_{p,r} [N/m^3]',true,climsr,zlims,rlims);

ax5=subplot(3,2,6);
plotvalue(ax5,Fdrd,'-mnn_n<\sigma_d v>u_r','F_{d,r} [N/m^3]',true,climsr,zlims,rlims);


linkaxes([ax0 ax1 ax2 ax3 ax4 ax5],'xy')


f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';


M.savegraph(f,sprintf('%s/%sforce_balancer_%i_art',folder,name,floor(mean(it))),papsize)


%% Azimuthal forces
f=figure('Name','azimuthal');


ax0=subplot(3,2,1);
plotvalue(ax0,duzd,'mn\partial_tu_{\theta}','F_{a,\theta} [N/m^3]',true,climsthet,zlims,rlims);

ax1=subplot(3,2, 2);
plotvalue(ax1,ones(size(FEzd)),'qnE_\theta', 'F_{E,\theta} [N/m^3]',true,climsthet,zlims,rlims);

ax2=subplot(3,2,3);
plotvalue(ax2,FBthetd,'qn(u_zB_r-u_rB_z)','F_{B,\theta} [N/m^3]',true,climsthet,zlims,rlims);

ax3=subplot(3,2,4);
plotvalue(ax3,FIthetd,'-mn(\bf{u\bullet\nabla})u_\theta','F_{i,\theta} [N/m^3]',true,climsthet,zlims,rlims);

ax4=subplot(3,2,5);
plotvalue(ax4,FPthetd,'-(\nablaP)_\theta','F_{p,\theta} [N/m^3]',true,climsthet,zlims,rlims);

ax5=subplot(3,2,6);
plotvalue(ax5,Fdthetd,'-mnn_n<\sigma_d v>u_\theta','F_{d,\theta} [N/m^3]',true,climsthet,zlims,rlims);



linkaxes([ax0 ax1 ax2 ax3 ax4 ax5],'xy')

f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
f.PaperSize=papsize;

M.savegraph(f,sprintf('%s/%sforce_balancethet_%i_art',folder,name,floor(mean(it))),papsize)

%% Axial forces
f=figure('Name','axial');
%zlim=[Fmin Fmax];

%zlim=[-15 15];

ax0=subplot(3,2,1);
plotvalue(ax0,duzd,'mn\partial_tu_{z}','F_{a,z} [N/m^3]',true,climsz,zlims,rlims);

ax1=subplot(3,2,2);
plotvalue(ax1,FEzd,'qnE_z', 'F_{E,z} [N/m^3]',true,climsz,zlims,rlims);

ax2=subplot(3,2,3);
plotvalue(ax2,FBzd,'-qn(u_\thetaB_r)','F_{B,z} [N/m^3]',true,climsz,zlims,rlims);

ax3=subplot(3,2,4);
plotvalue(ax3,FIzd,'-mn(\bf{u\bullet\nabla})u_z','F_{i,z} [N/m^3]',true,climsz,zlims,rlims);

ax4=subplot(3,2,5);
plotvalue(ax4,FPzd,'-(\nablaP)_z','F_{p,z} [N/m^3]',true,climsz,zlims,rlims);

ax5=subplot(3,2,6);
plotvalue(ax5,Fdzd,'-mnn_n<\sigma_d v>u_z','F_{d,z} [N/m^3]',true,climsz,zlims,rlims);

linkaxes([ax0 ax1 ax2 ax3 ax4 ax5],'xy')

f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
f.PaperSize=papsize;

M.savegraph(f,sprintf('%s/%sforce_balancez_%i_art',folder,name,floor(mean(it))),papsize)


    function plotvalue(ax,matrix,titl, clab, rm1strow, clims, zlims, rlims)
        if nargin < 5
            rm1strow=false;
        end
        if nargin>5
            cmin=clims(1);
            cmax=clims(2);
        else
            cmin=-inf;
            cmax=inf;
        end
        matrix(M.geomweight(:,:,1)<=0)=NaN;
        if rm1strow
            if log
                h=contourf(ax,M.zgrid*1e3,M.rgrid(1:end)*1e3,abs(matrix(1:end,:)),30,'linecolor','none');
            else
                h=contourf(ax,M.zgrid*1e3,M.rgrid(1:end)*1e3,matrix(1:end,:),30,'linecolor','none');
            end
        else
            if log
                h=contourf(ax,M.zgrid*1e3,M.rgrid*1e3,abs(matrix),30,'linecolor','none');
            else
                h=contourf(ax,M.zgrid*1e3,M.rgrid*1e3,matrix,30,'linecolor','none');
            end
        end
        %set(h,'edgecolor','none');
        hold on;
        if log
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, abs(matrix), densitycontour(1,2:end)*1e3, densitycontour(2,2:end)*1e3);
        else
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, matrix, densitycontour(1,2:end)*1e3, densitycontour(2,2:end)*1e3);
        end
        %plot3(ax,densitycontour(1,2:end)*1e3,densitycontour(2,2:end)*1e3,zpos,'-.','linewidth',2,'color',contourcolor)
        contour(M.zgrid*1e3,M.rgrid*1e3,mean(N,3),[contourscale contourscale]*maxdens,'-.','linewidth',2,'color',contourcolor);

        contour(M.zgrid*1e3,M.rgrid*1e3,M.geomweight(:,:,1),[0 0],'--','linewidth',2,'color','red');
        
        %plot(ax,geometriccontour(1,2:end)*1e3,geometriccontour(2,2:end)*1e3,'--','linewidth',2,'color','red')
        xlabel(ax,'z [mm]')
        ylabel(ax,'r [mm]')
        
        if (nargin>6 && ~isempty(zlims))
            xlim(ax,zlims*1e3)
        else
            xlim(ax,[M.zgrid(1) M.zgrid(end)]*1e3)
        end
        if (nargin>7 && ~isempty(rlims))
            ylim(ax,rlims*1e3)
        else
            if M.conformgeom
                ylim(ax,[M.rgrid(1) M.rgrid(rgridmax)]*1e3)
            else
                ylim(ax,[M.rgrid(1) M.rgrid(end)]*1e3)
            end
        end
        
        colormap(ax,'jet')
        c = colorbar(ax);
        c.Label.String='[N/m^{3}]';
        caxis(ax,[cmin cmax])
        if norm
            caxis(ax,[Fminr Fmaxr]);
        end
        if log
            set(ax,'colorscale','log')
            caxis(ax,'auto')
        end
        title(ax,clab(1:end-8))
        grid(ax, 'on')
        view(ax,2)
%         pos=ax.Position;
%         pos(3)=0.98*pos(3);
%         pos(4)=0.95*pos(4);
%         ax.Position=pos;
        
    end
end

