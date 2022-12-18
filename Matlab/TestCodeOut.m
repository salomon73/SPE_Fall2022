%% To work from home on local machine (macbook) 
%cd /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Inputs/Test_ions3    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))
Ions = espic2dhdf5('stable_dt_12.h5');

%% To work from ppb110 
addpath /home/sguincha/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('ions_vert_slice_dt10.h5');

%% Display particles data
dispespicParts(Ions);

%% Display fields data
dispespicFields(Ions);

%% Numerical parameters %%
Ions_mass = 1.67262192000000e-27; %3.347e-27;
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

phiA = 0;      % change for potinn when available                                                                                                
phiB = -20000;

rA = Ions.r_a; % change with RA and RB
rB = Ions.r_b;

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
phiA  = 0;
phiB  = -20000;

f = @(b,x) b(1) .* log(b(2).*x) + b(3);                           % Log Fit With Y-Offset
B = fminsearch(@(b) norm(Etild - f(b,Rtild)), [(phiA-phiB)/log(rB/rA), 1./rA, phiA-phiA]); % Initial guess
    plot(Rtild, 1.e-3*f(B,Rtild), 'r-', 'linewidth', 2)
    hold on 
    plot(sort(R0),1e-3*(phiA-phiB)*log(sort(R0)./rA)./log(rB/rA), 'g-', 'linewidth', 2)
    legend(strcat('$dt = $', num2str(Ions.dt)),strcat('y = ', num2str(1e-3*B(1)), ...
                                                      '$\log($', num2str(B(2)), ...
                                                      '$\cdot R_0) + $', num2str(B(3))), ...
                                                      strcat('y = ', ...
                                                      '$\Delta \phi \cdot \log(\frac{r}{r_a})/ \log(r_b/r_a)$' ), ...
                                                      'Location','northwest','Interpreter','latex');
            
                                                  
                                                  
%% subplot for energy collected at electrode

figure 
    subplot(1,2,1)
        xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
        ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
        legend(strcat('$dt = $', num2str(Ions.dt)),strcat('y = ', num2str(1e-3*B(1)), ...
                                                          '$\log($', num2str(B(2)), ...
                                                          '$\cdot R_0) + $', num2str(B(3))), ...
                                                          strcat('y = ', ...
                                                          '$\Delta \phi \cdot \log(\frac{r}{r_a})/ \log(r_b/r_a)$' ), ...
                                                          'Location','northwest','Interpreter','latex');
                                                          set(legend,'FontSize',18);
                                                          set (gca, 'fontsize', 22)
subplot(1,2,2)
        xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
        ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
        legend(strcat('$dt = $', num2str(Ions.dt)),strcat('y = ', num2str(1e-3*B(1)), ...
                                                          '$\log($', num2str(B(2)), ...
                                                          '$\cdot R_0) + $', num2str(B(3))), ...
                                                          strcat('y = ', ...
                                                          '$\Delta \phi \cdot \log(\frac{r}{r_a})/ \log(r_b/r_a)$' ), ...
                                                          'Location','northwest','Interpreter','latex');
                                                          set(legend,'FontSize',18);
                                                          set (gca, 'fontsize', 22)

%% Energy collected at electrode semilog %%

figure
    semilogx(R0,1.e-3*Energy, 'ko')
    hold on 
    semilogx(sort(R0),1e-3*(phiA-phiB)*log(sort(R0)./rA)./log(rB/rA), 'g-', 'linewidth', 2)
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E_{el}$ [keV]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), ...
           strcat('y = ', '$\Delta \phi \cdot \log(\frac{r}{r_a})/ \log(r_b/r_a)$' ), ...
                  'Location','northwest','Interpreter','latex', ...
                  'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)


%% Kinetic energy components %%

outliersE  = ~isoutlier(Energy./ER);
outliersEt = ~isoutlier(Ethet./ER);
outliersEtr = ~isoutlier(Energy./(Ethet+ER));

figure
    plot(sort(R0(outliersE)), sort(Energy(outliersE)./ER(outliersE)), 'r-', 'linewidth', 2)
    hold on
    plot(sort(R0(outliersEt)), sort(Ethet(outliersEt)./ER(outliersEt)), 'b-', 'linewidth', 2)
    hold on 
    plot(sort(R0(outliersEtr)), sort(Energy(outliersEtr)./(ER(outliersEtr)+ Ethet(outliersEtr))), 'k-','linewidth', 2)
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$E/E^{R}$, $E^{\theta}/E^{R}$ []', 'Interpreter', 'Latex')
    legend('$E/E^{R}$', '$E^{\theta}/E^{R}$','$E/(E^{R}+E^{\theta})$', 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    



%% Energy loss in electrode and yield%%

out = Stainless();
E   = out.E;
Eloss = out.Eloss;
element = out.element;

dR0 = 1e-3;
EnergyRange = 1e-6*[min(Energy), max(Energy)];
IndicesFit  = find(E>EnergyRange(1) & E<=EnergyRange(2)+dR0);
RangeOfInterest = E(IndicesFit);


figure
    semilogx(E,Eloss, 'r+-', 'linewidth', 2)
    hold on
    plot(EnergyRange(1)*ones(1,141),linspace(0,2500,141), 'k--', 'linewidth', 2)
    hold on 
    plot(EnergyRange(2)*ones(1,141),linspace(0,2500,141), 'k--', 'linewidth', 2)
    hold on
    plot(E(IndicesFit), Eloss(IndicesFit), 'bo-');

DegFit = 4;
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

energy_coord = linspace(0.001,0.02,1e5); % energy coordinate for extrapolated yield
Yield_pol  =  LambdaExp * polyval(P,energy_coord); % yield on this r coordinate


figure 
    plot(E(IndicesFit), Yield, 'ko', 'Linewidth', 2);
    hold on
    plot(energy_coord, Yield_pol, 'r-', 'linewidth', 2)
    ylabel('$\gamma(E)$ []', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [MeV]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$\gamma$ for tabulated $dE/dx$ values', ...
            strcat('Extrapolated $\gamma$ for',element), ...
            'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)


%% Plot of energy as a function of initial position R0 and yield %%
r = linspace(rA,rB, 10000);
E_r = 1.602*1e-19/(1.602*1e-19) * (-(phiB-phiA)*log(r./rA)./log(rB./rA) + phiA);
yield = LambdaExp  * polyval(P,1e-6*E_r);

figure
subplot(1,2,1)
    plot(r, 1e-6*E_r, 'k-', 'linewidth', 2)
    ylabel('$E(r)$ [MeV]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$r$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)
subplot(1,2,2)
    plot(E_r, yield, 'r-', 'linewidth', 2)
    hold on 
    plot(1e6*E(IndicesFit), Yield, 'ko', 'Linewidth', 2);
    ylabel('$\gamma(E)$ ', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [MeV]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)
    

%% Plot of potential phi %%
radial_coord = linspace(0.001, 0.01, 1000000);
phiR = -(phiB-phiA)*log(radial_coord./rA)./log(rB./rA) + phiA;

figure
    plot(radial_coord,phiR, 'k-', 'linewidth', 2)
    xlabel('$R_{0}$ [m]', 'Interpreter', 'Latex')
    ylabel('$\phi$ [V]', 'Interpreter', 'Latex')
    legend(strcat('$dt = $', num2str(Ions.dt)), ...
           'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    hold on 

%% find energy values (for theoretical yield)
b(1) = (phiA-phiB)/log(rB/rA);
b(2) = 1/rA;
b(3) = 0.0;
f = @(b,x) b(1) .* log(b(2).*x) + b(3);
Er3 = f(b,0.003);
Er5 = f(b,0.005);
Er8 = f(b,0.008);

[val, Ind] = min(abs(0.5*Er3 - 1e6*energy_coord));
gam3 = 2*Yield_pol(Ind) 
[val, Ind] = min(abs(0.5*Er5 - 1e6*energy_coord));
gam5 = 2*Yield_pol(Ind)  
[val, Ind] = min(abs(0.5*Er8 - 1e6*energy_coord));
gam8 = 2*Yield_pol(Ind)  
 
%% Lower energy model %%

gamma = compute_yield_potential(Ions, 'H')

%% Full energy range yield %%

energy_coord = linspace(0.001,0.02,1e3); % energy coordinate for extrapolated yield
Yield_pol  =  LambdaExp * polyval(P,energy_coord); % yield on this r coordinate
Hydrogen   = compute_yield_potential(Ions, 'H');
Helium     = compute_yield_potential(Ions, 'He');
Neon       = compute_yield_potential(Ions, 'Ne');
npoints    = 1000;
figure
    plot(energy_coord, Yield_pol, 'r-', 'linewidth', 2)
    ylabel('$\gamma(E)$', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [MeV]', 'interpreter', 'latex', 'Fontsize', 22)
    hold on 
    plot(linspace(0, 1e-3, npoints), cat(2,Hydrogen*ones(1,npoints-1),min(Yield_pol)), 'k-', 'linewidth', 2)
    hold on 
    plot(linspace(0, 1e-3, npoints), cat(2,Neon*ones(1,npoints-1),min(Yield_pol)), 'b-', 'linewidth', 2)
    hold on 
    plot(linspace(0, 1e-3, npoints), cat(2,Helium*ones(1,npoints-1),min(Yield_pol)), 'm-', 'linewidth', 2)
    legend(strcat('Extrapolated $\gamma$ for',element), ...
           '$\gamma_{p}^{H}$', '$\gamma_{p}^{Ne}$', '$\gamma_{p}^{He}$', ...
            'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22);

%% Find best fit for dE/dx curve - piecewise %%

index   = find(E<5e-2);
energie = E(index);
dEdx    = Eloss(index);

figure
plot(energie, dEdx);
hold on 
plot(energie, polyval(P,energie), 'k--')
hold on 
plot(linspace(0,5e-2, 1000),polyval(P,linspace(0,5e-2, 1000)), 'r--' )


%% Energy loss in electrode and yield%%

out = Aluminum();
E   = out.E;
Eloss = out.Eloss;
element = out.element;

dR0 = 3.5e-2;
EnergyRange = 1e-6*[min(Energy), max(Energy)];
IndicesFit  = find(E>EnergyRange(1) & E<=EnergyRange(2)+dR0);
RangeOfInterest = E(IndicesFit);

IndicesFit1 = find(E>EnergyRange(1) & E <=0.01);
IndicesFit2  = find(E>=0.01 & E<=0.02);
IndicesFit3     = find(E>=0.02 & E <= 0.03);
IndicesFit4     = find(E>=0.03 & E<=0.05);

deg1 = 3;
deg2 = 3;
deg3 = 3;
deg4 = 3;

P1 = polyfit(E(IndicesFit1), Eloss(IndicesFit1), deg1);
y1 = polyval(P1,E(IndicesFit1));

P2 = polyfit(E(IndicesFit2), Eloss(IndicesFit2), deg2);
y2 = polyval(P2,E(IndicesFit2));

P3    = polyfit(E(IndicesFit3), Eloss(IndicesFit3), deg3);
y3    = polyval(P3,E(IndicesFit3));

P4    = polyfit(E(IndicesFit4), Eloss(IndicesFit4), deg4);
y4    = polyval(P4,E(IndicesFit4));

figure
    semilogx(E,Eloss, 'r+-', 'linewidth', 2)
    hold on
    plot(E(IndicesFit1),y1, 'k-', 'linewidth', 1);
    hold on 
    plot(E(IndicesFit2),y2, 'g-', 'linewidth', 1);
    hold on 
    plot(E(IndicesFit3),y3, 'b-', 'linewidth', 1);
    hold on 
    plot(E(IndicesFit4),y4, 'y-', 'linewidth', 1);
    hold on
    plot(EnergyRange(1)*ones(1,141),linspace(0,2500,141), 'w--', 'linewidth', 1)
    hold on 
    plot(E(IndicesFit1(1))*ones(1,141),linspace(0,2500,141), 'g--', 'linewidth', 1)
    hold on 
    plot(E(IndicesFit2(1))*ones(1,141),linspace(0,2500,141), 'k--', 'linewidth', 1)
    hold on 
    plot(E(IndicesFit3(1))*ones(1,141),linspace(0,2500,141), 'b--', 'linewidth', 1)
    hold on 
    plot(E(IndicesFit4(1))*ones(1,141),linspace(0,2500,141), 'k--', 'linewidth', 1)
    hold on 
    plot(E(IndicesFit4(end))*ones(1,141),linspace(0,2500,141), 'k--', 'linewidth', 1)
    ylabel('$\frac{dE}{dx}$ [Mev/cm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$E$ [Mev]', 'interpreter', 'latex', 'Fontsize', 22)
    legend('$dE/dx$', strcat('fitting $P(E)$: $n=$',num2str(deg1)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg2)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg3)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg4)),...
           'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)

%% 
LambdaExp = 1e-3;
yield1 = LambdaExp * polyval(P1,E(IndicesFit1));
yield2 = LambdaExp * polyval(P2,E(IndicesFit2));
yield3 = LambdaExp * polyval(P3,E(IndicesFit3));
yield4 = LambdaExp * polyval(P4,E(IndicesFit4));
allindex = cat(2, IndicesFit1, IndicesFit2, IndicesFit3, IndicesFit4);
yieldTabulated = LambdaExp * Eloss(allindex);
DeltaX = E(IndicesFit1(1))-0; %slope
DeltaYH = yieldTabulated(IndicesFit(1))-Hydrogen;
DeltaYHe = yieldTabulated(IndicesFit(1))-Helium;
DeltaYNe = yieldTabulated(IndicesFit(1))-Neon;
X = linspace(0,DeltaX,10);
hydrogen_yield = DeltaYH/DeltaX*X + Hydrogen;
helium_yield   = DeltaYHe/DeltaX*X + Helium;
neon_yield     = DeltaYNe/DeltaX*X + Neon;

figure
    plot(E(IndicesFit1), yield1, 'r-', 'linewidth', 2)
    hold on 
    plot(E(IndicesFit2), yield2, 'k-', 'linewidth', 2)
    hold on 
    plot(E(IndicesFit3), yield3, 'b-', 'linewidth', 2)
    hold on 
    plot(E(IndicesFit4), yield4, 'g-', 'linewidth', 2)
    hold on 
    plot(E(allindex), yieldTabulated, 'ko')
    hold on 
    plot(X,hydrogen_yield,'-', 'linewidth', 2)
    hold on 
    plot(X,helium_yield,'-', 'linewidth', 2)
    hold on 
    plot(X,neon_yield,'-', 'linewidth', 2)
    legend(           strcat('fitting $P(E)$: $n=$',num2str(deg1)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg2)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg3)),...
                      strcat('fitting $P(E)$: $n=$',num2str(deg4)),...
                      '$\gamma$',...
                      'H','He', 'Ne',...
           'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)

%% PPB110 %%
electrons = espic2dhdf5('resultrestart_1e-12.h5');
%%
dt_el = electrons.dt(end); %last particle specie : electrons 
tpart_el = electrons.species(3).tpart;
vperp_el = electrons.species(3).Vperp(:,:);
vpar_el  = electrons.species(3).Vpar(:,:);
R_el = electrons.species(3).R(1,:);


%% PLOT GYRATION MOTION %%
figure
    plot(tpart_el,1e3*R_el, 'b-', 'linewidth', 2);
    ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$t$ [s]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)

%%
PlotParticleTrajectory(electrons.species(3), 1,1:200)


%% test IIEE module fortran90
nbparts=100;
AddedElectrons = Ions.species(1);
Relec = AddedElectrons.R(1:nbparts,:);
Thetelec = AddedElectrons.THET(1:nbparts,:);
% NEED TO BE SORTED TAKING ACCOUNT FOR PARALLELISATION RESULTS
% with PARTINDEX


%%
seuils_E = zeros(1,4);
seuils_E(1) = E(IndicesFit1(1))
seuils_E(2) = E(IndicesFit2(1))
seuils_E(3) = E(IndicesFit3(1))
seuils_E(4) = E(IndicesFit4(1))

%% FIGURE ION INDUCED TREATMENT %%

% electrode domain
Rgrid = linspace(0,10,10);
Zgrid = linspace(0,1,2);

[R,Z] = meshgrid(Rgrid, Zgrid);
vect_comp_Z = [ 0.2 ];
vect_comp_R = [0.3];

figure
    hold on 
    plot(linspace(0,10,10), ones(1,10), 'k-', 'linewidth', 2)
    hold on 
    p = pcolor(R,Z,ones(2,10));
    set(p, 'FaceAlpha', 0.4)
    shading interp
    hold on 
    plot(2,1.04, 'bo', 'MarkerSize', 10,  'MarkerFaceColor', 'b')
    hold on 
    plot(2,0.67, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    hold on 
    plot(5,1.2, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
    hold on 
    plot(5,0.8, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    hold on 
    plot(8,1.1, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
    hold on 
    plot(8,0.7, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    hold on 
    q = quiver(5,1.2,  vect_comp_Z, vect_comp_R);
    set(q, 'Color', 'k')
    hold on 
    vect_comp_Z = [-0.2 0 0.5];
    vect_comp_R = [0.4583 0.3 0.4583];
    q =quiver(2*ones(1,3),1.04*ones(1,3),  vect_comp_Z, vect_comp_R);
    set(q, 'Color', 'k')
    hold on 
    vect_comp_Z = [-0.3 0.4];
    vect_comp_R = [0.2  0.4583];
    q =quiver(8*ones(1,2),1.1*ones(1,2),  vect_comp_Z, vect_comp_R);
    set(q, 'Color', 'k')
    ylabel('$R$', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$Z$', 'interpreter', 'latex', 'Fontsize', 22)
    legend('boundary', 'cathode', '$e_1$', '$i_1$', '$e_2$', '$i_2$', '$e_3$', '$i_3$', ...
            '$v_0$', '$v_0$','$v_0$', 'Interpreter', 'latex')
    set(legend,'FontSize',20 );
    set (gca, 'fontsize', 22)
    set(gca, 'xtick', [])
    set(gca, 'xticklabel', [])
    set(gca, 'ytick', [])
    set(gca, 'yticklabel', [])

%% Function definitions %%

function gamma = compute_yield_potential(Ions, Neutral)

    alpha  =  0.032; % alpha parameter from Hagstrum model
    beta   =  0.78;  % Beta parameter from Hagstrum model
    out    =  Stainless();
    EF     =  out.EFermi; 
    Charge =  Ions.qe/1.602176620e-19;
    d      =  (Charge + 3.7)*1e-10; % See Hasselkamp for distance model (Angstrum) 
    WF     =  1.602176620e-19*(Ions.potinn-Ions.potout)*Ions.phinorm/log(Ions.r_b/Ions.r_a)*log((Ions.r_a+d)/Ions.r_a)-EF;

switch Neutral 
    case 'H'
        gamma = 0.0; %Ei = 13.6 < 2*EF
    case 'He'
        Ei = 24.58738;
        gamma = alpha*(beta*Ei-2*abs(WF));
    case 'Ne'
        Ei = 21.5646;
        gamma = alpha*(beta*Ei-2*abs(WF));
    otherwise 
        error('Neutral gas should be Hydrogen, Helium or Neon');
end

% end of function
end
