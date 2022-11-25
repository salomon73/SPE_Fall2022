addpath( '/home/sguincha/espic2d/matlab/')
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))

filename = 'Test_H_Al.h5';
ions = espic2dhdf5(filename);
dt = ions.dt;
time = dt*linspace(0,double(ions.nrun),ions.nrun+1);
%% 
figure
    plot(time, ions.nbparts, 'linewidth', 2)
    hold on 
    plot(time, ions.species(1).nbparts, 'linewidth', 2)
    xlabel('t [s]', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    legend('$n_i$', '$n_e$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

