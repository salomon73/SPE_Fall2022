%%
addpath '/home/sguincha/espic2d/matlab/'
%electrons1 = espic2dhdf5('RunT_1/resultrestart_5e-12.h5');
%electrons2 =  espic2dhdf5('resultrestart_5e-12.h5');

electrons = espic2dhdf5('/scratch/sguincha/SPE_Fall2022/Results/resultrestart_5e-12.h5');

%%
Nspecies = size(electrons.species());
R = electrons.species(end).R;
R = R(:,:);
%% plot electrons trajectories - parts1%%
val = 1:2
PlotParticleTrajectory(electrons.species(4), 1 , 1:2200)

%% 
PlotParticleTrajectory(electrons.species(3),1,1:2020)

%%
dt_el1 = electrons1.dt(end); %last particle specie : electrons 
tpart_el1 = electrons1.species(end).tpart;
dt_el2 = electrons2.dt(end); %last particle specie : electrons 
tpart_el2 = electrons2.species(end).tpart;
% vperp_el = electrons.species(end).Vperp;
% vperp_el = vperp_el(:,:);
% vpar_el  = electrons.species(end).Vpar;
% vpar_el = vpar_el(:,:);
R_el1 = electrons1.species(end).R;
R_el1 = R_el1(:,:);
R_el2 = electrons2.species(end).R;
R_el2 = R_el2(:,:);

%%
figure
    plot(tpart_el1,1e3*R_el1(1,:), 'b-', 'linewidth', 1);
    hold on
    plot(tpart_el2,1e3*R_el2(1,:), 'r-', 'linewidth', 1 )
    ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$t$ [s]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)

%%
LocalMinInd = islocalmin(1e3*R_el(1,:));

%% PLOT GYRATION MOTION %%
figure
    plot(tpart_el,1e3*R_el(1,:), 'b-', 'linewidth', 2);
    hold on
    plot(tpart_el(LocalMinInd), 1e3*R_el(1,LocalMinInd), 'r*' )
    ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$t$ [s]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)