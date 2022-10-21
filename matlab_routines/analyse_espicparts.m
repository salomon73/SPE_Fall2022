file='teststablegauss_13_traject2.h5';
%close all hidden;
if (~exist('M','var'))
    M.file=file;
end
M=load_espic2d(file,M,'all');


tmin=2;
tmax=length(M.ekin);
figure
plot(M.t0d(tmin:tmax),M.ekin(tmin:tmax),'o-',...
    M.t0d(tmin:tmax),M.epot(tmin:tmax),'d-',...
    M.t0d(tmin:tmax),M.etot(tmin:tmax),'h-',...
    M.t0d(tmin:tmax),M.eerr(tmin:tmax),'x--')
legend('ekin', 'epot', 'etot','eerr')
xlabel('Time [s]')
ylabel('Energies [J]')
grid on

figure
semilogy(M.t0d(tmin:tmax),abs(M.eerr(tmin:tmax)/M.etot(2)*100),'h-')
xlabel('t [s]')
ylabel('E_{err} %')
xlim([M.t0d(tmin),M.t0d(tmax)])
grid on


% Rlarmor=(M.me*sqrt(M.VR.^2+M.VTHET.^2)/(M.qe*M.B0));
% Rmax=M.R;
% Rmin=M.R-Rlarmor;
% [Isoutsideuprow,Isoutsideupcol]=find(Rmax-M.rgrid(end)>0);
% [Isoutsidedownrow,Isoutsidedowncol]=find(M.rgrid(1)-Rmin>0);

%dispespicParts(M)


f=figure();
tstudied=2;
legtext=sprintf("t=%2.1f - %2.1f [ns]",M.tpart(end-tstudied)*1e9,M.tpart(end)*1e9);
subplot(1,2,1)
H=M.H(:,1);
h1=histogram(H,20,'BinLimits',[min(H(:)) max(H(:))],'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
hold on
H=mean(M.H(:,end-tstudied:end),2);
h1=histogram(H,'BinWidth',h1.BinWidth,'DisplayName',legtext);
ylabel('counts')
xlabel('H [J]')
legend

subplot(1,2,2)
P=M.P(:,1);
h2=histogram(P,20,'BinLimits',[min(P(:)) max(P(:))],'DisplayName',sprintf("t=%2.3d [ns]",M.tpart(1)*1e9));
hold on
P=mean(M.P(:,end-tstudied:end),2);
histogram(P,'BinWidth',h2.BinWidth,'DisplayName',legtext);
ylabel('counts')
xlabel('P [kg\cdotm^2\cdots^{-1}]')
xlim(h2.BinLimits)
%clear P
%clear H
legend
f.PaperOrientation='landscape';
[~, name, ~] = fileparts(M.file);
xlim([0.95*h2.BinLimits(1) 1.05*h2.BinLimits(2)])
print(f,sprintf('%sParts_HP',name),'-dpdf','-fillpage')





BrillouinRatio=2*M.omepe^2/M.omece^2