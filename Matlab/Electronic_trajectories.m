%% Process electron trajectories for electrons generated at electrodes %%

path = '/scratch/sguincha';
directory = '/SPE_HRES_RUN/ResultsScan3/'; % Change result folder name accordingly
directory2 = 'SPE_Fall2022/RunsElectronsSlices/Results/';
directory3 = '/scratch/sguincha/SPE_Fall2022/Results/';
%-----------------------------------------
addpath(genpath(strcat(path,directory3)));
addpath(genpath('/home/sguincha/espic2d/matlab/'))
format long 

%% Create output structure containing electron infos %%

electrons = espic2dhdf5('resultrestart_5e-12.h5');
species    = size(electrons.species);
nbSpecies  = species(2);
Rel        = electrons.species(end).R;
R          = Rel(:,:);
nparts     = length(R(:,1));
npartalloc = 10;
tpart      = electrons.species(end).tpart;
ngyrations = 5;
%% RUN test - isparticle leaving electrode ? %%
for i = 1:npartalloc
    
    LocalMinInd(i,:) = islocalmin(1e3*R(i,:)); % test if each time step has local min of elect traject
    
    if isempty(find(LocalMinInd(i,:), 1,'first'))
        isfarenough(i) = 1;
        index(i)       = NaN;
        break
    else 
        index(i)         = find(LocalMinInd(i,:), 1,'first');  % find first local min index
        diff(i)          = R(i,index(i))-R(i,1);
        isfarenough(i)   = (diff(i)>1e-8);
    end
end

%% PLOT R trajectory if farenough ==1 %%

for ii = 1:npartalloc
    switch isfarenough(i)
        case 1
            figure
                plot(1e6*tpart,1e3*R(ii,:), 'b-', 'linewidth', 2);
                hold on
                plot(1e6*tpart(index(ii)), 1e3*R(ii,index(ii)), 'r*' )
                ylabel('$r_e$ [mm]', 'interpreter', 'latex','Fontsize', 22)
                xlabel('$t$ [$\mu$s]', 'interpreter', 'latex', 'Fontsize', 22)
                set (gca, 'fontsize', 22)
        case 0 
            disp('Careful ! After 1 Larmor gyration, particle may not be far enough ')
            
    end 
    
end

%% Plot electrons trajectories in vacuum vessel %%

% This loop ensures that we do not try to plot over time steps longer than
% existence time of the electron in the vacuum vessel 
for j = 1:npartalloc
    for k= 2:length(tpart)
        if R(j,k) == R(j,k-1)
            timeStep(j) = k-1;
        break 
        else 
            timeStep(j)=k;
        end
    end
end

    PlotParticleTrajectory(electrons.species(nbSpecies), 1:npartalloc, 1:min(timeStep(:)))

