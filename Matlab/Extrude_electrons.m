%% Process electron trajectories for electrons generated at electrodes %%

    path = '/scratch/sguincha';

    addpath(genpath('/home/sguincha/espic2d/matlab/'))
    addpath(genpath('/home/sguincha/SPE_Fall2022/Matlab/'))
    format long 

    electrons   = espic2dhdf5('resultfast.h5');
    
    %% Particles characteristics %%\
    
    
    % Numerical particles characteristics %    
    RelVn       = electrons.species(end).R;
    RVn         = RelVn(:,:);
    npartsVn    = length(RVn(:,1));
    nrun        = length(tpart);
    ngyrations  = 5;
    partindicesVn   = electrons.species(end).partindex(:,:);
    
    electronsVn    = electrons.species(end);
    LastTimeStepVn = zeros(1,npartsVn); % pre allocation for speed
    isfarenoughVn  = zeros(1,npartsVn); % pre allocation for speed
    indexVn        = zeros(1,npartsVn); % pre allocation for speed
    diffVn         = zeros(1,npartsVn); % pre allocation for speed
    LocalMinIndVn  = zeros(1,nrun);     % pre allocation for speed
    
    % Physical constants 
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    format long;
    
    % Energy characteristics 
    lowerBound = 0.1;
    upperBound = 10;
    nPoints    = 20;
    nEnergy    = nPoints;
    E = linspace(lowerBound, upperBound, nPoints);

    % Input electrons characteristics 
    nElectrons = length(tel);
    nblocks = nEnergy * nElectrons; 
    nparts  = nEnergy * nElectrons;
    

    
    %% RUN test - isparticle leaving electrode ? %%
    tic
    % NORMAL COMPONENTS LOOP %
    for i = 1:npartsVn
    
        LocalMinIndVn(i,:) = islocalmin(1e3*RVn(i,:));    % test if each time step has local min of elect traject

        if isempty(find(LocalMinIndVn(i,:), 1,'first')) % if no local min has been found: check if particle is recaptured
            isfarenoughVn(i) = 0;
            indexVn(i)       = NaN;

        else 
            indexVn(i)         = find(LocalMinIndVn(i,:), 1,'first');  % find first local min index
            diffVn(i)          = RVn(i,indexVn(i))-RVn(i,1);
            isfarenoughVn(i)   = (diffVn(i)>1e-8);
        end
    end

    % Loops to ensure that we do not try to plot over time steps longer than
    % existence time of the electrons in the vacuum vessel 

    for j = 1:npartsVn 
        for k= 2:nrun
            if RVn(j,k) == RVn(j,k-1)
                LastTimeStepVn(j) = k-1;
            break % find last living time step of electron(j)
            else 
                LastTimeStepVn(j)=k;
            end
        end
    end
    toc
    
    
    
    
%% Particles trajectories processing for given energies %% 

    nElectrons  = 25;
    nPoints     = 10;

    EnergyPartsIndices = zeros(1,nPoints);
    PositionSameCompo  = zeros(nPoints,nElectrons);
    

    % Gives first index of particle for energy E(ii)
    for ii = 1:nPoints
       EnergyPartsIndices(ii) = (ii-1)*nbPartsperEnergy+1; 
    end

    
    
%==============================================================================================================
% BELOW : COMPONENTS NORMAL TO ELECTRODE ONLY
%==============================================================================================================
    
%% Scan Normal component %% 

    energyVal = 2; % must be between 1 and nPoints = length(E)
    posVal    = 4; % must be between 1 and nElectrons
    
    PositionSameEnergVn  = zeros(nPoints,nElectrons); % all particles positions for a given energy
    
    % All positions corresponding to same energy %
    for ii = 1:nPoints
        
        PositionSameEnergVn(ii,:) = (1+nElectrons*(ii-1)):(1+nElectrons*(ii-1)+nElectrons -1);
        
    end
    EnergiesForSamePos = PositionSameEnergVn'; % energies for a given position are given by columns
    
    
%% Plot all particles with a given energy (normal v0) %%
    disp(strcat('E = ', num2str(E(energyVal)), ' eV'));
    disp('Initial positions: all');
    PlotParticleTrajectory(electronsVn,PositionSameEnergVn(energyVal,:),1:nrun)

%% Plot trajectories for all energy values at given position %%
    disp(strcat('(R0,Z0) = (', num2str(PartInfoVn(2,posVal)),',',num2str(PartInfoVn(3,posVal)), ')'));
    disp(strcat('E in [', num2str(E(1)),',' ,num2str(E(end)),  '] eV'));
    disp('V0 = Vr eR')
    PlotParticleTrajectory(electronsVn,EnergiesForSamePos(posVal,:),1:nrun)
    
    
    