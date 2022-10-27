%% Process electron trajectories for electrons generated at electrodes %%

path = '/scratch/sguincha';
directory = '/SPE_HRES_RUN/ResultsScan3/'; % Change result folder name accordingly
directory2 = 'SPE_Fall2022/RunsElectronsSlices/Results/';
directory3 = '/scratch/sguincha/SPE_Fall2022/Results/';
%-----------------------------------------
addpath(genpath(strcat(path,directory3)));
addpath(genpath('/home/sguincha/espic2d/matlab/'))
addpath(genpath('/home/sguincha/SPE_Fall2022/Matlab/'))
format long 

%% Create output structure containing electron infos %%

electrons   = espic2dhdf5('resultrestart_5e-12.h5');
species     = size(electrons.species);
nbSpecies   = species(2);
Rel         = electrons.species(end).R;
R           = Rel(:,:);
nparts      = length(R(:,1));
npartalloc  = 10;
tpart       = electrons.species(end).tpart;
ngyrations  = 5;
partindices = electrons.species(3).partindex(:,:);
velocity_struct = Temperature();
E            = velocity_struct.E;
LastTimeStep = zeros(1,npartalloc); % pre allocation for speed
isfarenough  = zeros(1,npartalloc); % pre allocation for speed
index        = zeros(1,npartalloc); % pre allocation for speed
diff         = zeros(1,npartalloc); % pre allocation for speed

%% RUN test - isparticle leaving electrode ? %%
for i = 1:npartalloc
    
    LocalMinInd(i,:) = islocalmin(1e3*R(i,:));    % test if each time step has local min of elect traject
    
    if isempty(find(LocalMinInd(i,:), 1,'first')) % if no local min has been found: check if particle is recaptured
        isfarenough(i) = 0;
        index(i)       = NaN;
        
    else 
        index(i)         = find(LocalMinInd(i,:), 1,'first');  % find first local min index
        diff(i)          = R(i,index(i))-R(i,1);
        isfarenough(i)   = (diff(i)>1e-8);
    end
end

% This loop (below) ensures that we do not try to plot over time steps longer than
% existence time of the electron in the vacuum vessel 

for j = 1:npartalloc 
    for k= 2:length(tpart)
        if R(j,k) == R(j,k-1)
            LastTimeStep(j) = k-1;
        break % find last living time step of electron(j)
        else 
            LastTimeStep(j)=k;
        end
    end
end

%% PLOT R trajectory if farenough ==1 %%

for ii = 1:npartalloc
    switch isfarenough(ii)
        case 1
            figure
                plot(1e6*tpart(1:LastTimeStep(ii)),1e3*R(ii,1:LastTimeStep(ii)), 'b-', 'linewidth', 2);
                hold on
                plot(1e6*tpart(index(ii)), 1e3*R(ii,index(ii)), 'r*' )
                ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
                xlabel('$t$ [$\mu$s]', 'interpreter', 'latex', 'Fontsize', 22)
                legend(strcat('E=',num2str(E(ii)), ' [eV]'), 'Location','northwest','Interpreter','latex');
                set (gca, 'fontsize', 22)
        case 0 
            disp(strcat('Careful ! After 1 Larmor gyration, particle may not be far enough for particle with E=', num2str(E(ii)), ' eV'))
            captured_part_indices(ii) = ii;
    end 
    
end
RecollectParts = nonzeros(captured_part_indices);

figure
    for kk = 1:length(RecollectParts)
        plot(1e6*tpart(1:LastTimeStep(RecollectParts(kk))+1),1e3*R(RecollectParts(kk),1:LastTimeStep(RecollectParts(kk))+1), '--', 'linewidth',0.5 );
        hold on
    end
    ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
    xlabel('$t$ [$\mu$s]', 'interpreter', 'latex', 'Fontsize', 22)
    legendstring = 'E=' + string(E(RecollectParts)) + ' [eV]';
    legend(legendstring, 'Location','northwest','Interpreter','latex');
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)
%% Plot R trajectory for recaptured particles %% 

figure 
    for ii = 3:npartalloc-1
        hold on 
        plot(1e6*tpart(1:100), 1e3*R(ii,1:100), 'linewidth', 1)
        ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
        xlabel('$t$ [$\mu$s]', 'interpreter', 'latex', 'Fontsize', 22)
        legendstring = 'E=' + string(E(3:npartalloc-1)) + ' [eV]';
        legend(legendstring, 'Location','northwest','Interpreter','latex');
        set(legend,'FontSize',18);
        set(gca,   'FontSize',22)
    end


%% Plot electrons trajectories in vacuum vessel %%

GoodParts = npartalloc - (find(isnan(index), 1, 'first')+1);

    PlotParticleTrajectory(electrons.species(nbSpecies), 1:GoodParts, 1:min(LastTimeStep(1:GoodParts)))
    PlotParticleTrajectory(electrons.species(nbSpecies), RecollectParts, 1:min(LastTimeStep(RecollectParts)))

