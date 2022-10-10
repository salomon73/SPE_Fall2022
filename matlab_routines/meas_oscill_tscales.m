
% t0=M.t0d;
% [t0,it]=unique(t0);
% npart=M.npart(it);
it=1:40:length(M.t2d);
t0=M.t2d(it);
N=M.N(:,:,it);
geomw=M.geomweight(:,:,1);
geomw(geomw<0)=0;
geomw(geomw>0)=1;

N=N.*geomw;
npart=squeeze(max(max(N,[],1),[],2));
[t0,it]=unique(t0);
npart=npart(it);


f=figure;
semilogy(t0,npart,'linewidth',1.5,'displayname','raw')
grid
hold on

npart=lowpass(npart,1/3e-7,1/(t0(2)-t0(1)));

[peaksmax,locsmax]=findpeaks(npart,t0,'minpeakdistance',1.2e-7);
[peaksmin,locsmin]=findpeaks(-npart,t0,'minpeakdistance',1.2e-7);

semilogy(t0,npart,'linewidth',1.5,'displayname','filtered')
grid
hold on
plot(locsmin,-peaksmin,'v','linewidth',1.5,'markersize',11)
plot(locsmax,peaksmax,'^','linewidth',1.5,'markersize',11)

peaksmin=-peaksmin;
peaksmax=peaksmax(2:end);
locsmax=locsmax(2:end);
if(locsmax(end)>locsmin(end))
    locsmax=locsmax(1:end-1);
    peaksmax=peaksmax(1:end-1);
end
if (peaksmin(1)<mean(peaksmin(2:end))/10)
peaksmin=peaksmin(2:end);
locsmin=locsmin(2:end);
end
peaksmin=peaksmin(1:end-1);
locsmin=locsmin(1:end-1);
if length(locsmin)>length(locsmax)
    nui=log(peaksmax./peaksmin(1:end-1))./(locsmax-locsmin(1:end-1));
    nul=log(peaksmax./peaksmin(2:end))./(locsmax-locsmin(2:end));
elseif length(locsmin)==length(locsmax)
    nui=log(peaksmax./peaksmin)./(locsmax-locsmin);
    nul=log(peaksmax(1:end-1)./peaksmin(2:end))./(locsmax(1:end-1)-locsmin(2:end));
else
    nui=log(peaksmax(1:end-1)./peaksmin)./(locsmax(1:end-1)-locsmin);
    nul=log(peaksmax(1:end-2)./peaksmin(2:end))./(locsmax(1:end-2)-locsmin(2:end));
end


nuimean=mean(nui)
nulmean=mean(nul)

%nuimean=1.8e8;
%nulmean=-2.6e8;

for i=1:length(peaksmax)
    t=[locsmax(i), locsmax(i)+2e-7];
    plot(t,peaksmax(i)*exp((t-t(1)).*nulmean),'--','linewidth',1.5)
end
for i=1:length(peaksmin)
    t=[locsmin(i), locsmin(i)+2e-7];
    plot(t,peaksmin(i)*exp((t-t(1)).*nuimean),'--','linewidth',1.5)
end

xlabel('time [ns]')
ylabel('N_e')
title(sprintf('\\Delta\\phi=%2.1fkV',(M.potout-M.potinn)*M.phinorm/1e3))
M.savegraph(f,sprintf('%s/%s_nuoscillsn',M.folder,M.name))