%%
if(isfile(sprintf('%s_DensPot_fft.mat',[M.folder,'/',M.name])))
    S=whos('-file',sprintf('%s_DensPot_fft.mat',[M.folder,'/',M.name]));
    if( ~exist('potfft','var') && sum(contains({S.name},'potfft')))
        load(sprintf('%s_DensPot_fft.mat',[M.folder,'/',M.name]),'potfft')
    end
    if( ~exist('densfft','var') && sum(contains({S.name},'densfft')))
        load(sprintf('%s_DensPot_fft.mat',[M.folder,'/',M.name]),'densfft')
    end
else
    densfft=M.N(:,:,:);
    potfft=M.pot(:,:,:);
    save(sprintf('%s_DensPot_fft',[M.folder,'/',M.name]),'densfft','potfft','-v7.3')
end
%%
fig=figure;
ax1=gca;
contourf(ax1,M.zgrid,M.rgrid,densfft(:,:,end));
hold on
[r,z]=find(densfft(:,:,end)~=0);
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= 'n [m^{-3}]';
view(ax1,2)
%set(ax1,'colorscale','log')

stepmin=100;
stepend=1500;%size(densfft,3);
papsize=[16 14];

%%
[x,y]=ginput(1);
nz=find(x>M.zgrid,1,'last');
nr=find(y>M.rgrid,1,'last');
nr=28;nz=143;

% decoupled
%nr=18; nz=127;
%nr=18; nz=76;

%nr=35; nz=192;
%stable
% nr=18; nz=63;
%nr=18; nz=38;
Z=(M.zgrid(nz)+M.zgrid(nz+1))/2;
R=(M.rgrid(nr)+M.rgrid(nr+1))/2;
Z=(M.zgrid(nz));
R=(M.rgrid(nr));
plot(Z,R,'rx','Markersize',12);

%decoupled
% nr=20;nz=129;
%nr=18; nz=127;
%nr=18; nz=76;
%nr=19; nz=127;
%nr=16;nz=128;

% stable
%nr=18; nz=63;
% nr=18; nz=38;

fig.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
fig.PaperUnits='centimeters';
fig.PaperSize=papsize;
print(fig,sprintf('%s_fftpoints',name),'-dpdf','-fillpage')
savefig(fig,sprintf('%s_fftpoints',name))

%% Compute standard deviation

figure;
ax1=gca;
surfc(ax1,M.zgrid,M.rgrid,std(densfft(:,:,stepmin:stepend),0,3));
hold on
[r,z]=find(densfft(:,:,end)~=0);
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= '\sigma_n [m^{-3}]';
%set(ax1,'colorscale','log')
view(ax1,2)

figure;
ax1=gca;
surfc(ax1,M.zgrid,M.rgrid,std(potfft(:,:,stepmin:stepend),0,3));
hold on
xlim(ax1,[M.zgrid(min(z)) M.zgrid(max(z))])
ylim(ax1,[M.rgrid(min(r)) M.rgrid(max(r))])
xlabel(ax1,'Z [m]')
ylabel(ax1,'R [m]')
c = colorbar(ax1);
c.Label.String= '\sigma_\phi [V]';
%set(ax1,'colorscale','log')
view(ax1,2)



%% Start the fourrier transform
dt=M.t2d(2)-M.t2d(1);

FS=1/dt;
L=stepend-stepmin;
f=FS*double(0:L/2)/L/1e9;

S=potfft(nr,nz,stepmin:stepend);

S=detrend(S(:),0);
S=S/max(S);
S=lowpass(S(:),.499*FS,FS);


P2=abs(fft(S(:)))/L;
P1=P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

tmin=stepmin*dt;
tmax=dt*size(densfft,3);

%% Compute the relevant frequencies
fce=M.qe/M.me/2/pi*M.Bz(nz,nr)/1e9;
N=densfft(:,:,end);
fpe=sqrt(mean(densfft(nr,nz,stepmin:end))*M.qe^2/(M.eps_0*M.me))/(2*pi*1e9);
fpemin=sqrt(min(min(N(N>0.01e14)))*M.qe^2/(M.eps_0*M.me))/(2*pi*1e9);
fpemax=sqrt(max(max(N))*M.qe^2/(M.eps_0*M.me))/(2*pi*1e9);
freplus=-sign(M.qsim)*fce/2*(1+sqrt(1-2*fpe^2/fce^2));
freminus=-sign(M.qsim)*fce/2*(1-sqrt(1-2*fpe^2/fce^2));
zIndices=mean(M.Presstens(6,nr,2:end-1,stepend-300:stepend),4)>0;
zIndices=find(zIndices);
lengths=Blocklength(mean(densfft(nr,:,stepmin:stepend),3)>1e10);
%zlength=M.width*2;
zlength=max(lengths)*M.dz;
v2=squeeze((M.Presstens(1,nr,zIndices,stepend-300:stepend)+M.Presstens(6,nr,zIndices,stepend-300:stepend))/M.me);
invdensv2=1./squeeze(densfft(nr,zIndices,stepend-300:stepend));
invdensv2(isinf(invdensv2))=0;
fz=sqrt(mean(mean(v2.*invdensv2)))/zlength/1e9;

%% Plot the fft
fig=figure('Name',sprintf('%s fft phi',M.name));
papsize=[14 16];

subplot(3,1,1)
pl=plot(f(1:end),P1(1:end),'linewidth',1.2,'Displayname','Spectrum');
title(sprintf('r=%3.2e[m],z=%3.2e[m], t=[%3.2e,%3.2e][s]',R,Z,tmin,tmax))
xlabel('f [GHz]','Interpreter','latex','fontsize',14)
ylabel('$\left|\hat{\Phi}(f)\right|$ [a.u.]','Interpreter','latex','fontsize',14)
grid on
ylimits=ylim();
hold on
plot(freminus*[1 1],ylimits,'k-.','Displayname','f_{re}^-','linewidth',1.2);
plot(freplus*[1 1],ylimits,'k-','Displayname','f_{re}^+','linewidth',1.2);
plot(fce*[1 1],ylimits,'r--','Displayname','f_{ce}','linewidth',1.2);
plot(fz*[1 1],ylimits,'k:','Displayname','f_{z}','linewidth',1.2);
plot(fpe*[1 1],ylimits,'r:','Displayname','local f_{pe}','linewidth',1.2);
fill([fpemin fpemax fpemax fpemin],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)],'r','FaceColor',[1 0 0]*0.5,'FaceAlpha',.3,'Displayname','f_{pe} range');
legend('location','eastoutside')
uistack(pl,'top');
ylim(ylimits)
xlim([min(f),1.2*max(fpemax)])
%%
subplot(3,1,2)
pl=plot(f(1:end),P1(1:end),'linewidth',1.2,'Displayname','Spectrum');
xlabel('f [GHz]','Interpreter','latex','fontsize',14)
ylabel('$\left|\hat{\Phi}(f)\right|$ [a.u.]','Interpreter','latex','fontsize',14)
grid on
%xlim([0 10])
ylimits=ylim();
hold on
plot(freminus*[1 1],ylimits,'k-.','Displayname','f_{re}^-','linewidth',1.2);
plot(freplus*[1 1],ylimits,'k-','Displayname','f_{re}^+','linewidth',1.2);
plot(fce*[1 1],ylimits,'r--','Displayname','f_{ce}','linewidth',1.2);
plot(fz*[1 1],ylimits,'k:','Displayname','f_{z}','linewidth',1.2);
plot(fpe*[1 1],ylimits,'r:','Displayname','f_{pe}','linewidth',1.2);
fill([fpemin fpemax fpemax fpemin],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)],'r','FaceColor',[1 0 0]*0.5,'FaceAlpha',.3,'Displayname','f_{pe}');
uistack(pl,'top');
ylim(ylimits)
xlim([min(f),max(f)])


subplot(3,1,3)
Signal=potfft(nr,nz,1:stepend);
%Signal=M.N(nr,nz,stepmin:end);
pl=plot(dt*((1:length(Signal))-1),Signal(:));
xlabel('t [s]','Interpreter','latex','fontsize',14)
ylabel('$\Phi(t)$ [V]','Interpreter','latex','fontsize',14)
ylimits=ylim();
%ylimits=[-56 -54];
hold on
x=[stepmin-1 length(Signal)-1 length(Signal)-1 stepmin-1]*dt;
y=[ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
fill(x,y,'r','FaceColor',[1 1 1]*0.5,'FaceAlpha',.3)
ylim(ylimits)
uistack(pl,'top');
%fig.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
fig.PaperUnits='centimeters';
fig.PaperSize=papsize;

print(fig,sprintf('%s_fftPhi_R%dZ%d',name,nr,nz),'-dpdf','-fillpage')
savefig(fig,sprintf('%s_fftPhi_R%dZ%d',name,nr,nz))


S=densfft(nr,nz,stepmin:stepend);
S=detrend(S(:),0);
S=S/max(S);
S=lowpass(S(:),.499*FS,FS);
P2=abs(fft(S(:)))/L;
P1=P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

fig=figure('Name',sprintf('%s fft dens',M.name));
%%
subplot(3,1,1)
pl=plot(f(1:end),P1(1:end),'linewidth',1.2,'Displayname','Spectrum');
title(sprintf('r=%3.2e[m],z=%3.2e[m], t=[%3.2e,%3.2e][s]',R,Z,tmin,tmax))
xlabel('f [GHz]','Interpreter','latex','fontsize',14)
ylabel('$\left|\hat{n}(f)\right|$ [a.u.]','Interpreter','latex','fontsize',14)
grid on
%xlim([0 10])
ylimits=ylim();
hold on
plot(freminus*[1 1],ylimits,'k-.','Displayname','f_{re}^-','linewidth',1.2);
plot(freplus*[1 1],ylimits,'k-','Displayname','f_{re}^+','linewidth',1.2);
plot(fce*[1 1],ylimits,'r--','Displayname','f_{ce}','linewidth',1.2);
plot(fz*[1 1],ylimits,'k:','Displayname','f_{z}','linewidth',1.2);
plot(fpe*[1 1],ylimits,'r:','Displayname','local f_{pe}','linewidth',1.2);
fill([fpemin fpemax fpemax fpemin],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)],'r','FaceColor',[1 0 0]*0.5,'FaceAlpha',.3,'Displayname','range f_{pe}');
legend('location','eastoutside')
uistack(pl,'top');
xlim([min(f),1.2*max(fpemax)])
%%
subplot(3,1,2)
pl=plot(f(1:end),P1(1:end),'linewidth',1.2,'Displayname','Spectrum');
xlabel('f [GHz]','Interpreter','latex','fontsize',14)
ylabel('$\left|\hat{n}(f)\right|$ [a.u.]','Interpreter','latex','fontsize',14)
grid on
%xlim([0 10])
ylimits=ylim();
hold on
plot(freminus*[1 1],ylimits,'k-.','Displayname','f_{re}^-','linewidth',1.2);
plot(freplus*[1 1],ylimits,'k-','Displayname','f_{re}^+','linewidth',1.2);
plot(fce*[1 1],ylimits,'r--','Displayname','f_{ce}','linewidth',1.2);
plot(fz*[1 1],ylimits,'k:','Displayname','f_{z}','linewidth',1.2);
plot(fpe*[1 1],ylimits,'r:','Displayname','f_{pe}','linewidth',1.2);
fill([fpemin fpemax fpemax fpemin],[ylimits(1) ylimits(1) ylimits(2) ylimits(2)],'r','FaceColor',[1 0 0]*0.5,'FaceAlpha',.3,'Displayname','f_{pe}');
uistack(pl,'top');
xlim([min(f),max(f)])
%%
subplot(3,1,3)
Signal=densfft(nr,nz,1:stepend);
%Signal=M.N(nr,nz,stepmin:end);
pl=plot(dt*((1:length(Signal))-1),Signal(:));
xlabel('t [s]','Interpreter','latex','fontsize',14)
ylabel('$n(t)$ [m$^{-3}$]','Interpreter','latex','fontsize',14)
ylimits=ylim();
hold on
x=[stepmin-1 length(Signal)-1 length(Signal)-1 stepmin-1]*dt;
y=[ylimits(1) ylimits(1) ylimits(2) ylimits(2)];
fill(x,y,'r','FaceColor',[1 1 1]*0.5,'FaceAlpha',.3)
uistack(pl,'top');
%fig.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
fig.PaperUnits='centimeters';
fig.PaperSize=papsize;
ylim(ylimits)
print(fig,sprintf('%s_fftdens_R%dZ%d',name,nr,nz),'-dpdf','-fillpage')
savefig(fig,sprintf('%s_fftdens_R%dZ%d',name,nr,nz))