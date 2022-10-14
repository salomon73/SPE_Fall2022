%% To work from home on local machine (macbook) 
cd /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Inputs/Test_ions3
addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('stable_dt_11.h5');

%% To work from ppb110 
cd /home/sguincha/SPE_Fall2022/Inputs/Test_ions3
addpath /home/sguincha/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('stable_dt_11.h5');
%% Display particles data
dispespicParts(Ions);

%% Display fields data
dispespicFields(Ions);

%% Follow the particles

Ions_mass = 3.347e-27;
Vr = Ions.VR;
Vz = Ions.VZ;
Vt = Ions.VTHET;
VR = Vr(:,:);
VZ = Vz(:,:);
VT = Vt(:,:);

R = Ions.R;
Z = Ions.Z;

Ra = 0.001;
Rb = 0.01;
Za = -0.32;
Z_b = 0.32;

dt = Ions.dt;
nt = R.nt;
t  = 5*dt*linspace(0,200,nt);
dr = 1.75e-3;

npart = R.nparts;
mat = zeros(npart,nt);
R = R(:,:);
Z = Z(:,:);

[R0 I] = sort(R(Ions.partindex(:,1),1), 'descend');

POS0 = [I, R0];

%%
tic 
parfor ii = 1:npart
    posR{ii}   = R(Ions.partindex(:,:)==ii);
    index(ii)  = length(posR{ii});
    mask       = Ions.partindex(:,index(ii))==ii;
    Energy(ii) = (1/1.602e-19)*0.5*Ions_mass*(VR(mask,index(ii))^2+VZ(mask,index(ii))^2+VT(mask,index(ii))^2);
    Ethet(ii)  = (1/1.602e-19)*0.5*Ions_mass*VT(mask,index(ii))^2;
    ER(ii)     = (1/1.602e-19)*0.5*Ions_mass*VR(mask,index(ii))^2;
    R0(ii)     = posR{ii}(1);
end
toc



%%

figure
    plot(R0,1.e-3*Energy, 'ko')
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    hold on 

[Rsort,sortId] = sort(R0);
Esort = Energy(sortId);
Etild = Esort(1:3:end);
Rtild = Rsort(1:3:end);

f = @(b,x) b(1) .* log(b(2).*x) + b(3);                % Exponential Fit With Y-Offset
B = fminsearch(@(b) norm(Etild - f(b,Rtild)), [1e4, 0.5, 73200]);         % Estimate Parameters

    plot(Rtild, 1.e-3*f(B,Rtild), 'r-', 'linewidth', 2)
    legend(strcat('$dt = $', num2str(Ions.dt)),strcat('y = ', num2str(1e-3*B(1)), ...
                                                      '$\log($', num2str(B(2)), ...
                                                      '$\cdot R_0) + $', num2str(B(3))), ...
                                                      'Location','northwest','Interpreter','latex');
            
%%

figure
    semilogx(R0,1.e-3*Energy, 'ko')
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)


%%
outliersE  = ~isoutlier(Energy./ER);
outliersEt = ~isoutlier(Ethet./ER);
outliersEtr = ~isoutlier(Energy./(Ethet+ER));

figure
    plot(R0(outliersE), Energy(outliersE)./ER(outliersE), 'ro')
    hold on
    plot(R0(outliersEt), Ethet(outliersEt)./ER(outliersEt), 'bo')
    hold on 
    plot(R0(outliersEtr), Energy(outliersEtr)./(ER(outliersEtr)+ Ethet(outliersEtr)), 'ko')
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E/E^{R}$, $E^{\theta}/E^{R}$ []', 'Interpreter', 'Latex')
    legend('$E/E^{R}$', '$E^{\theta}/E^{R}$','$E/(E^{R}+E^{\theta})$', 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    



%%
E0 = [.001 .002 .003 .004 .005 .007 .009 .01];
E1 = [.02 .03 .04 .05 .06 .07 .08 .09 ];
E2 = [.1 .2 .3 .4 .5 .6 .7 .8 .9 ];
E3 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5];
E4 = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95];
E5 = [100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 10000];

    E = cat(2, E0, E1, E2, E3, E4, E5);


Eloss0 = [210.09 296.26 362.38 418.12 467.23 564.17 647.61 685.32];
Eloss1 = [958.52 1103.4 1177.2 1209.7 1218.5 1213.4 1200.5 1183.1];
Eloss2 = [1163.3 969.92 838.72 747.47 678.39 622.57 575.31 530.43 493.76];
Eloss3 = [463.10 359.80 297.61 255.35 224.62 201.17 182.63 167.55 155.02 ...
           144.44 135.40 127.52 120.60 114.46 108.99 104.07 99.622 95.579];
Eloss4 = [91.884 67.098 53.524 44.910 38.903 34.461 31.033 28.303 26.075 ...
           24.219 22.649 21.301 20.132 19.106 18.202 17.396 16.673 16.021 ];
Eloss5 = [15.430 11.584 9.5581 8.3666 7.5440 6.9545 6.5131 6.1718 5.9011 ... 
           5.6822 5.5024 5.3527 5.2267 5.1197 5.0281 4.9492 4.8809 4.8124 4.7694];

    Eloss = cat(2, Eloss0, Eloss1, Eloss2, Eloss3, Eloss4, Eloss5);

dR0 = 1e-2;
EnergyRange = 1e-6*[min(Energy), max(Energy)];
IndicesFit  = find(E>EnergyRange(1) & E<=EnergyRange(2)+dR0);
RangeOfInterest = E(IndicesFit);


figure
    semilogx(E,Eloss, 'r+-', 'linewidth', 2)
    hold on 
    plot(EnergyRange(1)*ones(1,141),linspace(0,1400,141), 'k--', 'linewidth', 2)
    hold on 
    plot(EnergyRange(2)*ones(1,141),linspace(0,1400,141), 'k--', 'linewidth', 2)
    hold on
    plot(E(IndicesFit), Eloss(IndicesFit), 'bo-');

DegFit = 3;
P = polyfit(E(IndicesFit), Eloss(IndicesFit),DegFit);
y = polyval(P,E(IndicesFit));

    hold on 
    plot(E(IndicesFit),y, 'k-');
    ylabel('$\frac{dE}{dx}$ [Mev/cm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [Mev]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$dE/dx$', '$E(R_0^{min})$','$E(R_0^{max})$','fitted values',strcat('fitting $P(E)$: $n=$',num2str(DegFit)), 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)

LambdaExp = 1.e-3; % cm/MeV;
Yield = LambdaExp * y;

figure 
    plot(E(IndicesFit), Yield, 'ro-', 'Linewidth', 2);
    ylabel('$\gamma(E(R_0))$ []', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$R_0$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$\gamma$ for tabulated $dE/dx$ values', 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)











%% TRASH %%

%% Collection at electrode using dr interval to model contact %%
Ions_mass = 3.347e-27;
Vr = Ions.VR;
Vz = Ions.VZ;
Vt = Ions.VTHET;
VR = Vr(:,:);
VZ = Vz(:,:);
VT = Vt(:,:);

R = Ions.R;
Z = Ions.Z;

Ra = 0.001;
Rb = 0.01;
Za = -0.32;
Z_b = 0.32;

dt = Ions.dt;
nt = R.nt;
t  = 5*dt*linspace(0,200,nt);
dr = 1.75e-3;

npart = R.nparts;
mat = zeros(npart,nt);
R = R(:,:);

for jj = 1:npart
    for ii = 1:nt
        if R(2,ii) < Ra+dr
            mat(jj,ii) = ii; % Trouve l'indice des particles et du temps auquel elles se trouvent pres de/sur l'électrode 
        else 
            mat(jj,ii) = 0;
            
        end
    end
end

indices = zeros(1,2);
for jj =1:npart
    for ii = 1:nt
        t_ind = find(mat(jj,:)); % trouve le premier indice de temps où il y'a contact pour la particule jj
        indices(jj,1) = jj;
        indices(jj,2) = t_ind(1);
    end
end
k=150; % Index of particle of interest
EnergyIon_n = (1/1.602e-19)*0.5*Ions_mass*(VR(k,indices(k,2)).^2+VZ(k,indices(k,2)).^2+VT(k,indices(k,2)).^2);

for n=1:npart
    E(n) = (1/1.602e-19)*0.5*Ions_mass*(VR(n,indices(n,2)).^2+VZ(n,indices(n,2)).^2+VT(n,indices(n,2)).^2);
end

figure
plot(linspace(1,1000,1000),E, 'k+')


%% TO correct
for ii = 1:length(R0)                               % go along all particles array
    for jj = 1:nt                                   % go along all time steps
        if Ions.nbparts(jj) ~= 0                    % check if there are still particles in the cell
            if ismember(ii,Ions.partindex(:,jj))    % check if particle ii still in cell
                continue
            else 
                index(ii)  = jj;                 
                %mask       = Ions.partindex(:,index(ii))==jj;
                %Energy(ii) = (1/1.602e-19)*0.5*Ions_mass*(VR(mask,index(ii)).^2+VZ(mask,index(ii)).^2+VT(mask,index(ii)).^2);
                %EThet(ii)  = (1/1.602e-19)*0.5*Ions_mass*VT(ii,jj)^2;
                %ER(ii)  = (1/1.602e-19)*0.5*Ions_mass*VR(ii,jj)^2;
                break 
            end
        end
    end
end

