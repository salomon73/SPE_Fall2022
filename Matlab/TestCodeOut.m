%% To work from home on local machine (macbook) 
cd /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Inputs/Test_ions3
addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))
Ions = espic2dhdf5('stable_dt_11.h5');

%% To work from ppb110 
cd /home/sguincha/SPE_Fall2022/Inputs/Test_ions3
addpath /home/sguincha/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('stable_dt_11.h5');
%% Display particles data
dispespicParts(Ions);

%% Display fields data
dispespicFields(Ions);

%% Numerical parameters %%

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

%% Compute energy components and initial radial distance for ions %%

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



%% Energy collected at electrode with log fit %%

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

f = @(b,x) b(1) .* log(b(2).*x) + b(3);                           % Log Fit With Y-Offset
B = fminsearch(@(b) norm(Etild - f(b,Rtild)), [1e4, 0.5, 73200]); % Initial guess

    plot(Rtild, 1.e-3*f(B,Rtild), 'r-', 'linewidth', 2)
    hold on 
    plot(sort(R0),1e-3*(phiA-phiB)*log(sort(R0)./rA)./log(rB/rA), 'g-', 'linewidth', 2)
    legend(strcat('$dt = $', num2str(Ions.dt)),strcat('y = ', num2str(1e-3*B(1)), ...
                                                      '$\log($', num2str(B(2)), ...
                                                      '$\cdot R_0) + $', num2str(B(3))), ...
                                                      strcat('y = ', ...
                                                      '$\Delta \phi \cdot \log(\frac{r}{r_a})/ \log(r_b/r_a)$' ), ...
                                                      'Location','northwest','Interpreter','latex');
            
%% Energy collected at electrode semilog %%

figure
    semilogx(R0,1.e-3*Energy, 'ko')
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)


%% Kinetic energy components %%

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
    



%% Energy loss in electrode and yield%%

out = Aluminum();
E   = out.E;
Eloss = out.Eloss;
element = out.element;

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

DegFit = 6;
P = polyfit(E(IndicesFit), Eloss(IndicesFit),DegFit);
y = polyval(P,E(IndicesFit));

    hold on 
    plot(E(IndicesFit),y, 'k-');
    ylabel('$\frac{dE}{dx}$ [Mev/cm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [Mev]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$dE/dx$', '$E(R_0^{min})$','$E(R_0^{max})$','fitted values',...
            strcat('fitting $P(E)$: $n=$',num2str(DegFit)), ...
            'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)

LambdaExp = 1.e-3; % cm/MeV;
Yield = LambdaExp * y; 

rcoord = linspace(0.001,0.02,1e3); % r coordinate for extrapolated yield
Yield_pol  =  LambdaExp * polyval(P,rcoord); % yield on this r coordinate


figure 
    plot(E(IndicesFit), Yield, 'ko', 'Linewidth', 2);
    hold on
    plot(rcoord, Yield_pol, 'r-', 'linewidth', 2)
    ylabel('$\gamma(E(R_0))$ []', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$R_0$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$\gamma$ for tabulated $dE/dx$ values', ...
            strcat('Extrapolated $\gamma$ for',element), ...
            'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)




%% Plot of potential phi %%
phiA = 0;
phiB = -20000;

rA = Ions.r_a;
rB = Ions.r_b;
phiR = -(phiB-phiA)*log(R0./rA)./log(rB./rA) + phiA;

figure
    plot(R0,phiR, 'ko')
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$\phi$ [V]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), ...
           'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    hold on 

%% Lower energy model %%












%% TRASH %%

for i =1:10 % for loop to contract section in one line
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

end

