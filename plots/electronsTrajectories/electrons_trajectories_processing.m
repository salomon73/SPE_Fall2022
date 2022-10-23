addpath '/home/sguincha/espic2d/matlab/'
electrons = espic2dhdf5('results2parts.h5');

%% plot electrons trajectories - parts1%%

PlotParticleTrajectory(electrons.species(3), 1:2, 1:100)

%% 
PlotParticleTrajectory(electrons.species(3),2,1:100)

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
