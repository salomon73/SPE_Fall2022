%% generates paths and output structure %% 
addpath( '/home/sguincha/espic2d/matlab/')
addpath(genpath('/Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Matlab/Data/'))

filename = 'Test_H_SS.h5';
ions = espic2dhdf5(filename);
%%
dt = ions.dt;
time = dt*linspace(0,double(ions.nrun),ions.nrun+1);

%% plot number of particles over time %% 
figure
    plot(ions.species(2).tpart, ions.species(2).nbparts, 'linewidth', 2)
    hold on 
    plot(ions.species(1).tpart, ions.species(1).nbparts, 'linewidth', 2)
    hold on 
    plot(ions.tpart, ions.nbparts, 'linewidth', 2)
    xlabel('t [s]', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    legend('$n_i$', '$n_e$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

figure
    plot(ions.species(2).nbparts, 'linewidth', 2)
    hold on 
    plot(ions.species(1).nbparts, 'linewidth', 2)
    hold on 
    plot(ions.nbparts, 'linewidth', 2)
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
    plot(Z(2,76:end), 'k-', 'linewidth', 2)
    xlabel('nsteps', 'Interpreter', 'Latex') 
    ylabel('Z', 'Interpreter', 'Latex')
    set (gca, 'fontsize', 22)
    

%%
me = 9.1e-31;
vE = sqrt((2/me)*2*1.602e-19)

% %% 
% elec_id = ions.species(4).partindex;
% % elec_id = elec_id(:100,:);

Relec = ions.species(4).R;
Relec = Relec(1:2000,:);
Zelec = ions.species(4).Z;
Zelec = Zelec(1:2000,:);
nrun = 20000;
%%
threshold = 1e-6;
lostparts = zeros(1,10);
%% 
for ii = 600 : 610
   for jj = 2:nrun 
      if ((Relec(ii:jj) ==  Relec(ii,jj-1)) &  (Relec(ii, jj-1) ~= 0))
          lostparts(ii-599) = 0; % lost
      break
      else 
          lostparts(ii - 599) = 1; % not lost
      end 
   end
end


%% figure

figure
plot(Zelec(1500,2000:20000) , Relec(1500, 2000:20000), 'k*')
    xlabel('Z', 'Interpreter', 'Latex') 
    ylabel('R', 'Interpreter', 'Latex')
    set (gca, 'fontsize', 22)

    
    
 %%
 elec_id = ions.species(4).partindex;
 elec_id = elec_id(1:2000,:);
 
 
 
 %% Plot particles 
 PlotParticleTrajectory(ions.species(4), 2:3, 18000:20000) % for extrude_iiee_peak_density
 % save the above trajectory and superimpose with the potential well 
 ions.display2Dpotentialwell(0)
 
 
%% Test Code IIEE - Horizontal slices statistics %%

ionsSS = espic2dhdf5('Test_H_SS.h5');
ionsCu = espic2dhdf5('Test_H_Cu.h5');
ionsAl = espic2dhdf5('Test_H_Al.h5');
% colors
% '#77AC30' '#D95319' '#0072BD'
%%
figure
    p1 = plot(ionsSS.nbparts, 'linewidth', 2);
    set(p1(1), 'Color','k')
    hold on 
    p2 = plot(ionsSS.species(1).nbparts, 'linewidth', 2);
    set(p2(1), 'Color', '#D95319')
    hold on 
    p3 = plot(ionsCu.species(1).nbparts, 'linewidth', 2);
    set(p3(1), 'Color', '#77AC30')
    hold on 
    p4 = plot(ionsAl.species(1).nbparts, 'linewidth', 2);
    set(p4(1), 'Color', '#0072BD')
    xlabel('nsteps', 'Interpreter', 'Latex') 
    ylabel('nparts', 'Interpreter', 'Latex')
    xlim([0 400])
    legend('$n_i$', '$n_e^{SS}$' ,'$n_e^{Cu}$' ,'$n_e^{Al}$' ,'Location','best','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
    axes('position',[.55 .55 .10 .30])
    box on % put box around new pair of axes
    indexOfInterest = 280:320; % range of t near perturbation   
    p2 = plot(ionsSS.species(1).nbparts(indexOfInterest), 'linewidth', 2);
    set(p2(1), 'Color', '#D95319')
    hold on 
    p3 = plot(ionsCu.species(1).nbparts(indexOfInterest), 'linewidth', 2);
    set(p3(1), 'Color', '#77AC30')
    hold on 
    p4 = plot(ionsAl.species(1).nbparts(indexOfInterest), 'linewidth', 2);
    set(p4(1), 'Color', '#0072BD')
    axis tight
    set (gca, 'fontsize', 22)
    xticks([290 300])
    xticklabels({'290' '300'})
