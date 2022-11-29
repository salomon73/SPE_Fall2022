%% generates paths and output structure %% 
addpath( '/home/sguincha/espic2d/matlab/')
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))

filename = 'Test_V0_H_SS.h5';
ions = espic2dhdf5(filename);
dt = ions.dt;
time = dt*linspace(0,double(ions.nrun),ions.nrun+1);

%% plot number of particles over time %% 
figure
    plot(time, ions.nbparts, 'linewidth', 2)
    hold on 
    plot(time, ions.species(1).nbparts, 'linewidth', 2)
    xlabel('t [s]', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    legend('$n_i$', '$n_e$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

figure
    plot( ions.nbparts, 'linewidth', 2)
    hold on 
    plot(ions.species(1).nbparts, 'linewidth', 2)
    xlabel('nsteps', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    legend('$n_i$', '$n_e$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

%% Collect information about generated electrons %%
R = ions.species(1).R;
R = R(:,:);
R = R(1:25,:);
Z = ions.species(1).Z;
Z = Z(:,:);
Z = Z(1:25,:);

%% Checks that initial electron velocity is normal %% 
figure
    plot(Z(2,76:end), R(2,76:end), 'k-', 'linewidth', 2)
    xlabel('Z', 'Interpreter', 'Latex') 
    ylabel('R', 'Interpreter', 'Latex')
    set (gca, 'fontsize', 22)
   
figure
    plot(Z(2,76:end), 'k-')
