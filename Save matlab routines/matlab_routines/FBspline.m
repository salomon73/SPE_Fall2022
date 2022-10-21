function [Forces, Continuity]=FBspline(M,Forces, Continuity, Density,norm,log,parper,contourscale, zlims, rlims, clims)
%ForceBalance Show the radial force balance
%   Plot the three radial forces for the given time-step it or averaged
%   over the range of time steps defined in it

it=Forces.it;
switch nargin
    case 4
        parper=false;
        contourscale=0.6;
        norm=false;
        log=false;
        zlims=[-inf inf];
    case 5
        parper=false;
        log=false;
        contourscale=0.6;
        zlims=[-inf inf];
    case 6
        parper=false;
        contourscale=0.6;
        zlims=[-inf inf];
    case 7
        contourscale=0.6;
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



%% Radial forces
f=figure('Name','radial');


ax0=subplot(4,2,1);
plotvalue(ax0,durd,'Mass flux time derivative', 'mn\partialu_{r}/\partialt [Nm^{-3}]',true,climsr,zlims,rlims);

ax1=subplot(4,2,2);
plotvalue(ax1,FErd,'Electric force density', 'F_{E,r} [Nm^{-3}]',true,climsr,zlims,rlims);

ax2=subplot(4,2,3);
plotvalue(ax2,FBrd,'Magnetic force density','F_{B,r} [Nm^{-3}]',true,climsr,zlims,rlims);

ax3=subplot(4,2,4);
plotvalue(ax3,FIrd,'Inertial force density','F_{i,r} [Nm^{-3}]',true,climsr,zlims,rlims);

ax4=subplot(4,2,5);
plotvalue(ax4,FPrd,'Pressure force density','F_{p,r} [Nm^{-3}]',true,climsr,zlims,rlims);

ax5=subplot(4,2,6);
plotvalue(ax5,Fdrd,'Drag force density','F_{d,r} [Nm^{-3}]',true,climsr,zlims,rlims);

ax6=subplot(4,2,[7,8]);
plotvalue(ax6,totforcerd,'Tot force density','F_{tot,r} [Nm^{-3}]',true,climsr,zlims,rlims);


linkaxes([ax0 ax1 ax2 ax3 ax4 ax5 ax6],'xy')
if parper
    sgtitle(f,sprintf('Perp forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
else
sgtitle(f,sprintf('Radial forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
end
f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%sforce_balancer_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%sforce_balancer_%i',name,floor(mean(it))))

%% Radial relative forces
f=figure('Name','radial relative');

color_limits=[-1 1]*2;
referencer=FErd;

ax1=subplot(3,2,1);
plotvalue(ax1,FErd./referencer,'Relative Electric force density', 'F_{Er}/F_{Er}',true,color_limits,zlims,rlims);

%zlim=[Fmin Fmax];
zlim=[-inf inf];

ax2=subplot(3,2,2);
plotvalue(ax2,FBrd./referencer,'Relative Magnetic force density','F_{Br}/F_{Er}',true,color_limits,zlims,rlims);

ax3=subplot(3,2,3);
plotvalue(ax3,FIrd./referencer,'Relative Inertial force density','F_i/F_{Er}',true,color_limits,zlims,rlims);

ax4=subplot(3,2,4);
plotvalue(ax4,FPrd./referencer,'Relative Pressure force density','F_{pr}/F_{Er}',true,color_limits,zlims,rlims);

ax5=subplot(3,2,5);
plotvalue(ax5,Fdrd./referencer,'Relative Drag force density','F_{dr}/F_{Er}',true,color_limits,zlims,rlims);

ax6=subplot(3,2,6);
plotvalue(ax6,totforcerd./referencer,'Total force density','F_{totr}/F_{Er}',true,color_limits,zlims,rlims);


linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')
if parper
    sgtitle(sprintf('Relative perp forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
else
sgtitle(sprintf('Relative radial forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
end
f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%s_relatforce_balancer_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%s_relatforce_balancer_%i',name,floor(mean(it))))


%% Azimuthal relative forces
f=figure('Name','azimuthal relative');
color_limits=[-1 1]*2;
referencethet=FBthetd;

ax1=subplot(3,2,1);
plotvalue(ax1,duthetd./referencethet,'Mass flux time derivartive','mn\partialu_\theta/\partialt/F_{B,\theta}',true,color_limits,zlims,rlims);

ax2=subplot(3,2,2);
plotvalue(ax2,FBthetd./referencethet,'Magnetic force density','F_{B,\theta}/F_{B,\theta}',true,color_limits,zlims,rlims);

ax3=subplot(3,2,3);
plotvalue(ax3,FIthetd./referencethet,'Inertial force density','F_{i,\theta}/F_{B,\theta}',true,color_limits,zlims,rlims);

ax4=subplot(3,2,4);
plotvalue(ax4,FPthetd./referencethet,'Pressure force density','F_{p,\theta}/F_{B,\theta}',true,color_limits,zlims,rlims);

ax5=subplot(3,2,5);
plotvalue(ax5,Fdthetd./referencethet,'Drag force density','F_{d,\theta}/F_{B,\theta}',true,color_limits,zlims,rlims);

ax6=subplot(3,2,6);
plotvalue(ax6,totforcethetd./referencethet,'Total force density','F_{tot,\theta}/F_{B,\theta}',true,color_limits,zlims,rlims);



linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')

sgtitle(sprintf('Relative azimuthal forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%srelatforce_balancethet_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%srelatforce_balancethet_%i',name,floor(mean(it))))


%% Azimuthal forces
f=figure('Name','azimuthal');



ax1=subplot(3,2,1);
plotvalue(ax1,duthetd,'Mass flux time derivartive','mn\partialu_\theta/\partialt [Nm^{-3}]',true,climsthet,zlims,rlims);

ax2=subplot(3,2,2);
plotvalue(ax2,FBthetd,'Magnetic force density','F_{B,\theta} [Nm^{-3}]',true,climsthet,zlims,rlims);

ax3=subplot(3,2,3);
plotvalue(ax3,FIthetd,'Inertial force density','F_{i,\theta} [Nm^{-3}]',true,climsthet,zlims,rlims);

ax4=subplot(3,2,4);
plotvalue(ax4,FPthetd,'Pressure force density','F_{p,\theta} [Nm^{-3}]',true,climsthet,zlims,rlims);

ax5=subplot(3,2,5);
plotvalue(ax5,Fdthetd,'Drag force density','F_{d,\theta} [Nm^{-3}]',true,climsthet,zlims,rlims);

ax6=subplot(3,2,6);
plotvalue(ax6,totforcethetd,'Total force density','F_{tot,\theta} [Nm^{-3}]',true,climsthet,zlims,rlims);


linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')

sgtitle(sprintf('Azimuthal forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%sforce_balancethet_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%sforce_balancethet_%i',name,floor(mean(it))))

%% Axial forces
f=figure('Name','axial');
%zlim=[Fmin Fmax];

%zlim=[-15 15];

ax0=subplot(4,2,1);
plotvalue(ax0,duzd,'Mass flux time derivartive','mn\partialu_z/\partialt [Nm^{-3}]',true,climsz,zlims,rlims);

ax1=subplot(4,2,2);
plotvalue(ax1,FEzd,'Electric force density', 'F_{E,z} [Nm^{-3}]',true,climsz,zlims,rlims);

ax2=subplot(4,2,3);
plotvalue(ax2,FBzd,'Magnetic force density','F_{B,z} [Nm^{-3}]',true,climsz,zlims,rlims);

ax3=subplot(4,2,4);
plotvalue(ax3,FIzd,'Inertial force density','F_{i,z} [Nm^{-3}]',true,climsz,zlims,rlims);

ax4=subplot(4,2,5);
plotvalue(ax4,FPzd,'Pressure force density','F_{p,z} [Nm^{-3}]',true,climsz,zlims,rlims);

ax5=subplot(4,2,6);
plotvalue(ax5,Fdzd,'Drag force density','F_{d,z} [Nm^{-3}]',true,climsz,zlims,rlims);

ax6=subplot(4,2,[7,8]);
plotvalue(ax6,totforcezd,'Tot force density','F_{tot,z} [Nm^{-3}]',true,climsz,zlims,rlims);

linkaxes([ax0 ax1 ax2 ax3 ax4 ax5 ax6],'xy')

if parper 
    sgtitle(sprintf('Parallel forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
else
sgtitle(sprintf('Axial forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
end

f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%sforce_balancez_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%sforce_balancez_%i',name,floor(mean(it))))

%% Axial relative forces
f=figure('Name','axial relative');
ax1=subplot(3,2,1);
color_limits=[-1 1]*2;
referencez=FEzd;
plotvalue(ax1,FEzd./referencez,'Relative Electric force density', 'F_{Ez}/F_{Ez}',true,color_limits,zlims,rlims);

%zlim=[Fmin Fmax];
zlim=[-inf inf];

ax2=subplot(3,2,2);
plotvalue(ax2,FBzd./referencez,'Relative Magnetic force density','F_{Bz}/F_{Ez}',true,color_limits,zlims,rlims);

ax3=subplot(3,2,3);
plotvalue(ax3,FIzd./referencez,'Relative Inertial force density','F_{Iz}/F_{Ez}',true,color_limits,zlims,rlims);

ax4=subplot(3,2,4);
plotvalue(ax4,FPzd./referencez,'Relative Pressure force density','F_{Pz}/F_{Ez}',true,color_limits,zlims,rlims);

ax5=subplot(3,2,5);
plotvalue(ax5,Fdzd./referencez,'Relative Drag force density','F_{dz}/F_{Ez}',true,color_limits,zlims,rlims);

ax6=subplot(3,2,6);
plotvalue(ax6,totforcezd./referencez,'Total force density','F_{totz}/F_{Ez}',true,color_limits,zlims,rlims);


linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'xy')

if parper
    sgtitle(sprintf('Relative Parallel forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
else
sgtitle(sprintf('Relative axial forces t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
end
f.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f.PaperUnits='centimeters';
papsize=[20 29];
f.PaperSize=papsize;
print(f,sprintf('%s_relatforce_balancez_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%s_relatforce_balancez_%i',name,floor(mean(it))))

%% Forces sums
f2=figure('Name','total');
ax4=subplot(4,1,1);
maxforcerd=1;%cat(3,FErd,FBrd,FPrd,FIrd);
%maxforcerd=max(maxforcerd,[],3);
if parper
plotvalue(ax4,totforcerd./maxforcerd,'Total Perp force density','F [Nm^{-3}]',true,climsr,zlims,rlims);
else
    plotvalue(ax4,totforcerd./maxforcerd,'Total radial force density','F [Nm^{-3}]',true,climsr,zlims,rlims);
end

ax5=subplot(4,1,2);
maxforcethetd=1;%cat(3,FEzd,FBzd,FPzd,FIzd);
%maxforcezd=max(maxforcezd,[],3);
%zlims=[-15 15];
plotvalue(ax5,totforcethetd./maxforcethetd,'Total azimuthal force density','F [Nm^{-3}]',true,climsthet,zlims,rlims);


ax6=subplot(4,1,3);
maxforcezd=1;%cat(3,FEzd,FBzd,FPzd,FIzd);
%maxforcezd=max(maxforcezd,[],3);
%zlims=[-15 15];
if parper
    plotvalue(ax6,totforcezd./maxforcezd,'Total par force density','F [Nm^{-3}]',true,climsz,zlims,rlims);
else
    plotvalue(ax6,totforcezd./maxforcezd,'Total axial force density','F [Nm^{-3}]',true,climsz,zlims,rlims);
end
ax7=subplot(4,1,4);
surface(ax7,M.zgrid*1e3,M.rgrid*1e3,mean(N,3),'edgecolor','none');
hold on;
zlevel=interp2(M.zgrid*1e3, M.rgrid*1e3, mean(N,3), densitycontour(1,2:end)*1e3, densitycontour(2,2:end)*1e3);
plot3(ax7,densitycontour(1,2:end)*1e3,densitycontour(2,2:end)*1e3,zlevel,'-.','linewidth',2,'color',contourcolor)
leng=size(geometriccontour,2);
strt=1;
while leng>1
    span=geometriccontour(2,strt);
    plot(ax7,geometriccontour(1,strt+1:strt+span)*1e3,geometriccontour(2,strt+1:strt+span)*1e3,'--','linewidth',2,'color','red')
    leng=leng-span-1;
    strt=strt+span+1;
end


colormap(ax7,'parula')
xlabel(ax7,'z [m]')
ylabel(ax7,'r [m]')
title(ax7,'Density')
c = colorbar(ax7);
c.Label.String= 'n[m^{-3}]';
view(ax7,2)
linkaxes([ax4 ax5 ax6 ax7])

if (nargin>8 && ~isempty(zlims))
    xlim(ax7,zlims*1e3)
else
    xlim(ax7,[M.zgrid(1) M.zgrid(end)]*1e3)
end
if (nargin>9 && ~isempty(rlims))
    ylim(ax7,rlims*1e3)
else
    if M.conformgeom
        ylim(ax7,[M.rgrid(1) M.rgrid(rgridmax)]*1e3)
    else
        ylim(ax7,[M.rgrid(1) M.rgrid(end)]*1e3)
    end
end

sgtitle(sprintf('Forces sum t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
f2.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f2.PaperUnits='centimeters';
papsize=[20 29];
f2.PaperSize=papsize;
print(f2,sprintf('%sforce_dens_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f2,sprintf('%sforce_dens_%i',name,floor(mean(it))))


%% Continuity
if( ~isempty(Continuity))
f3=figure('Name','Continuity');
it=Continuity.it;
t=M.t2d(it);

ax7=subplot(4,1,1);
plotvalue(ax7,mean(Continuity.dndt,3),'\partialn/\partialt','\partialn/\partialt [m^{-3}s^{-1}]',true)


zlim=[max([min(min(mean(Continuity.ndivu(2:end,:,:),3))),min(min(mean(Continuity.ugradn(2:end,:,:),3)))]) min([max(max(mean(Continuity.ndivu(2:end,:,:),3))),max(max(mean(Continuity.ugradn(2:end,:,:),3)))])];

ax8=subplot(4,1,2);
plotvalue(ax8,mean(sum(Continuity.ndivu,4),3),'n\nabla(u)','n\nabla(u)[m^{-3}s^{-1}]',true,zlim)

ax9=subplot(4,1,3);
plotvalue(ax9,mean(sum(Continuity.ugradn,4),3),'u\nabla(n)','u\nabla(n)[m^{-3}s^{-1}]',true,zlim)

ax10=subplot(4,1,4);
plotvalue(ax10,mean(sum(Continuity.ndivu,4)+sum(Continuity.ugradn,4)+Continuity.dndt,3),'Total','total[m^{-3}s^{-1}]',true,zlim)

sgtitle(sprintf('Continuity t=[%1.2g-%1.2g]s n_e=%1.2g m^{-3} Ns=%d',M.t2d(min(it)),M.t2d(max(it)),double(maxdens),length(it)))
f3.PaperOrientation='portrait';
[~, name, ~] = fileparts(M.file);
f3.PaperUnits='centimeters';
papsize=[20 19];
f3.PaperSize=papsize;
print(f3,sprintf('%scontinuity_%i',name,floor(mean(it))),'-dpdf','-fillpage')
savefig(f,sprintf('%scontinuity_%i',name,floor(mean(it))))
end

    function plotvalue(axis,matrix,titl, clab, rm1strow, clims, zlims, rlims)
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
                h=surface(axis,M.zgrid*1e3,M.rgrid(1:end)*1e3,abs(matrix(1:end,:)));
            else
                h=surface(axis,M.zgrid*1e3,M.rgrid(1:end)*1e3,matrix(1:end,:));
            end
        else
            if log
                h=surface(axis,M.zgrid*1e3,M.rgrid*1e3,abs(matrix));
            else
                h=surface(axis,M.zgrid*1e3,M.rgrid*1e3,matrix);
            end
        end
        set(h,'edgecolor','none');
        hold on;
        if log
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, abs(matrix), densitycontour(1,2:end)*1e3, densitycontour(2,2:end)*1e3);
        else
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, matrix, densitycontour(1,2:end)*1e3, densitycontour(2,2:end)*1e3);
        end
        plot3(axis,densitycontour(1,2:end)*1e3,densitycontour(2,2:end)*1e3,zpos,'-.','linewidth',2,'color',contourcolor)
        
        
        if log
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, abs(matrix), geometriccontour(1,2:end)*1e3, geometriccontour(2,2:end)*1e3);
        else
            zpos=interp2(M.zgrid*1e3, M.rgrid*1e3, matrix, geometriccontour(1,2:end)*1e3, geometriccontour(2,2:end)*1e3);
        end
        leng=size(geometriccontour,2);
        strt=1;
        while leng>1
            span=geometriccontour(2,strt);
            plot(axis,geometriccontour(1,strt+1:strt+span)*1e3,geometriccontour(2,strt+1:strt+span)*1e3,'--','linewidth',2,'color','red')
            leng=leng-span-1;
            strt=strt+span+1;
        end
        xlabel(axis,'z [mm]')
        ylabel(axis,'r [mm]')
        if (nargin>6 && ~isempty(zlims))
            xlim(axis,zlims*1e3)
        else
            xlim(axis,[M.zgrid(1) M.zgrid(end)]*1e3)
        end
        if (nargin>7 && ~isempty(rlims))
            ylim(axis,rlims*1e3)
        else
            if M.conformgeom
                ylim(axis,[M.rgrid(1) M.rgrid(rgridmax)]*1e3)
            else
                ylim(axis,[M.rgrid(1) M.rgrid(end)]*1e3)
            end
        end
        
        colormap(axis,'jet')
        c = colorbar(axis);
        c.Label.String=clab;
        caxis(axis,[cmin cmax])
        if norm
            caxis(axis,[Fminr Fmaxr]);
        end
        if log
            set(axis,'colorscale','log')
            caxis(axis,'auto')
        end
        title(axis,titl)
        grid(axis, 'on')
        view(axis,2)
    end
end

