%% Define the constants %%
kB = (25.7/298)*1e-3; % eV/K

%% Energy range for T range %%
E = linspace(0.1, 10, 5);
T = 1/kB*E;
%%
figure
    plot(tpart_el1,1e3*R_el1(1,:), 'b-', 'linewidth', 1);
    hold on
    plot(tpart_el2,1e3*R_el2(1,:), 'r-', 'linewidth', 1 )
    ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$t$ [s]', 'interpreter', 'latex', 'Fontsize', 22)
    set (gca, 'fontsize', 22)