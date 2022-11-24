addpath( '/home/sguincha/espic2d/matlab/')
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))

filename = 'Test_He.h5';
ions = espic2dhdf5(filename);

%% 
figure
    plot(ions.species(3).nbparts, 'linewidth', 2)
    hold on 
    plot(ions.species(4).nbparts, 'linewidth', 2)
    xlabel('nsteps', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    legend('$n_i$', '$n_e$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

