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

    electrons   = espic2dhdf5('resultrestart_5e-12.h5');

%%

    % Electrons characteristics % 
    % -> electrons initialised with normal velocity
    % -> electrons with varying components along (r,z)
    % Read PartInfoV0
    load('PartInfos_10_25.mat')
    
    
    % Added species characteristics % 
    species     = size(electrons.species);
    nbSpecies   = species(2);
    RelVn       = electrons.species(end-2).R;
    RelVneg     = electrons.species(end-1).R;
    RelV0       = electrons.species(end).R;
    RVn         = RelVn(:,:);
    RVneg       = RelVneg(:,:);
    RV0         = RelV0(:,:);
    npartsVn    = length(RVn(:,1));
    npartsVneg  = length(RVneg(:,1));
    npartsV0    = length(RV0(:,1));
    tpart       = electrons.species(end).tpart;
    nrun        = length(tpart);
    ngyrations  = 5;
    partindicesVn   = electrons.species(end-2).partindex(:,:);
    partindicesVneg = electrons.species(end-1).partindex(:,:);
    partindicesV0   = electrons.species(end).partindex(:,:);
    nbAddedSpecies  = nbSpecies-3;

    % Energy characteristics %
    lowerBound = 0.1;
    upperBound = 20;
    nPoints = 10;
    E   = linspace(lowerBound, upperBound, nPoints);
    
    
    % Particles characteristics % 
    electronsVn    = electrons.species(end-2);
    LastTimeStepVn = zeros(1,npartsVn); % pre allocation for speed
    isfarenoughVn  = zeros(1,npartsVn); % pre allocation for speed
    indexVn        = zeros(1,npartsVn); % pre allocation for speed
    diffVn         = zeros(1,npartsVn); % pre allocation for speed
    LocalMinIndVn  = zeros(1,nrun);     % pre allocation for speed
    
    electronsV0Neg   = electrons.species(end-1);
    LastTimeStepVNeg = zeros(1,npartsVneg); % pre allocation for speed
    isfarenoughVNeg  = zeros(1,npartsVneg); % pre allocation for speed
    indexVNeg        = zeros(1,npartsVneg); % pre allocation for speed
    diffVNeg         = zeros(1,npartsVneg); % pre allocation for speed
    LocalMinIndVNeg  = zeros(1,nrun);       % pre allocation for speed
    
    
    electronsV0    = electrons.species(end);
    LastTimeStepV0 = zeros(1,npartsV0); % pre allocation for speed
    isfarenoughV0  = zeros(1,npartsV0); % pre allocation for speed
    indexV0        = zeros(1,npartsV0); % pre allocation for speed
    diffV0         = zeros(1,npartsV0); % pre allocation for speed
    LocalMinIndV0  = zeros(1,nrun); % pre allocation for speed
    
    % Tests before entering run section %
    try 
        npartsVn   == length(PartInfoVn(1,:));
        npartsVneg == length(PartInfoVneg(1,:));
        npartsV0   == length(PartInfoV0(1,:));
 
    catch
        
        warning('Error: different number of particles from output file and from initialisation information')
        
    end
    
    
    
    
%% RUN test - isparticle leaving electrode ? %%
    tic
    % POSITIVE COMPONENTS LOOP %
    for i = 1:npartsV0
    
        LocalMinIndV0(i,:) = islocalmin(1e3*RV0(i,:));    % test if each time step has local min of elect traject

        if isempty(find(LocalMinIndV0(i,:), 1,'first')) % if no local min has been found: check if particle is recaptured
            isfarenoughV0(i) = 0;
            indexV0(i)       = NaN;

        else 
            indexV0(i)         = find(LocalMinIndV0(i,:), 1,'first');  % find first local min index
            diffV0(i)          = RV0(i,indexV0(i))-RV0(i,1);
            isfarenoughV0(i)   = (diffV0(i)>1e-8);
        end
    end
    % NEGATIVE COMPONENTS LOOP %
    for i = 1:npartsVneg
    
        LocalMinIndVNeg(i,:) = islocalmin(1e3*RVneg(i,:));    % test if each time step has local min of elect traject

        if isempty(find(LocalMinIndVNeg(i,:), 1,'first')) % if no local min has been found: check if particle is recaptured
            isfarenoughVNeg(i) = 0;
            indexVNeg(i)       = NaN;

        else 
            indexVNeg(i)         = find(LocalMinIndVNeg(i,:), 1,'first');  % find first local min index
            diffVNeg(i)          = RV0(i,indexVNeg(i))-RVneg(i,1);
            isfarenoughVNeg(i)   = (diffVNeg(i)>1e-8);
        end
    end
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

    for j = 1:npartsV0 
        for k= 2:nrun
            if RV0(j,k) == RV0(j,k-1)
                LastTimeStepV0(j) = k-1;
            break % find last living time step of electron(j)
            else 
                LastTimeStepV0(j)=k;
            end
        end
    end
    
    for j = 1:npartsVneg 
        for k= 2:nrun
            if RVneg(j,k) == RVneg(j,k-1)
                LastTimeStepVNeg(j) = k-1;
            break % find last living time step of electron(j)
            else 
                LastTimeStepVNeg(j)=k;
            end
        end
    end
    
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

    nComponents = 6;
    nElectrons  = 25;
    nPoints     = 10;

    nbPartsperEnergy   = nComponents*nElectrons;
    EnergyPartsIndices = zeros(1,nPoints);
    PositionSameCompo  = zeros(nPoints,nElectrons);
    PositionSameEnerg  = zeros(1,nComponents*nElectrons);
    
    PosAllCompoPerPart = zeros(nPoints,nComponents);

    % Gives first index of particle for energy E(ii)
    for ii = 1:nPoints
       EnergyPartsIndices(ii) = (ii-1)*nbPartsperEnergy+1; 
    end
    
%==============================================================================================================
% BELOW : INITIAL VELOCITY COMPONENTS VARYING FROM VR TO VZ PURELY
%==============================================================================================================

%% Plot particles trajectories for energy value given by energyVal for V0 %%

    energyVal = 5; % must be between 1 and nPoints = length(E) - energy value
    compoVal  = 1;  % must be between 1 and nComponents - (vr0,vz0) index
    posVal    = 2;  % must be between 1 and nElectrons  - electron initial position
    
    % Find all indices with same (vr,vz) for given value %
    for ii = 1:nPoints 
        
        PositionSameCompo(ii,:) = (1+nElectrons*(compoVal-1)) + (nElectrons*nComponents)*(ii-1) : (1+nElectrons*(compoVal-1) + nElectrons-1  ) + (nElectrons*nComponents)*(ii-1);
        
    end
    
    % find all indices for given energy value %
    for ii = 1:nElectrons*nComponents 
        
        PositionSameEnerg(ii) = ii + nPoints*(energyVal-1);
        
    end
    
    % find all components for same particle (psoition) %
    for ii =1:nPoints 
       
        for jj = 1:nComponents 
            
            PosAllCompoPerPart(ii,jj) = nElectrons*(jj-1)+posVal +nElectrons*nComponents*(ii-1); % eventually change nElec*nCompo = 72
            
        end 
    end
    
    [C, IEnerg, ICompo] = intersect(PositionSameEnerg,PositionSameCompo(energyVal,:));
  
    
%% Plot all particles for a given energy (by step of 3) %% 
    PlotParticleTrajectory(electronsV0,EnergyPartsIndices(energyVal):3:EnergyPartsIndices(energyVal+1)-1, 1:nrun)

%% Plot all particles with given initial components for fixed E %%
    disp('All particles with same initial velocity components for a given energy')
    disp(strcat('(V0R,V0Z) = (',num2str(PartInfoVneg(5,PositionSameCompo(energyVal,1))), ',' ...
                               ,num2str(PartInfoVneg(4,PositionSameCompo(energyVal,1))) , ')'))
    disp(strcat('E = ', num2str(E(energyVal)), ' eV'));
    PlotParticleTrajectory(electronsV0,PositionSameCompo(energyVal,:),1:nrun)
    
%% Plot all components for a given particle position %%
    PlotParticleTrajectory(electronsV0, PosAllCompoPerPart(energyVal,:), 1:nrun)
    
    

    
    
    
    
%==============================================================================================================
% BELOW : INITIAL VELOCITY COMPONENTS VARYING FROM VR TO VZ PURELY WITH
%         NEGATIVE VZ COMPONENTS 
%==============================================================================================================

%% Plot particles trajectories for energy value given by energyVal for V0 %%

    energyVal = 2; % must be between 1 and nPoints = length(E) - energy value
    compoVal  = 2;  % must be between 1 and nComponents - (vr0,vz0) index
    posVal    = 2;  % must be between 1 and nElectrons  - electron initial position
    
    % Find all indices with same (vr,vz) for given value %
    for ii = 1:nPoints 
        
        PositionSameCompo(ii,:) = (1+nElectrons*(compoVal-1)) + (nElectrons*nComponents)*(ii-1) : (1+nElectrons*(compoVal-1) + nElectrons-1  ) + (nElectrons*nComponents)*(ii-1);
        
    end
    
    % find all indices for given energy value %
    for ii = 1:nElectrons*nComponents 
        
        PositionSameEnerg(ii) = ii + nPoints*(energyVal-1);
        
    end
    
    % find all components for same particle (psoition) %
    for ii =1:nPoints 
       
        for jj = 1:nComponents 
            
            PosAllCompoPerPart(ii,jj) = nElectrons*(jj-1)+posVal +nElectrons*nComponents*(ii-1); % eventually change nElec*nCompo = 72
            
        end 
    end
    
    [C, IEnerg, ICompo] = intersect(PositionSameEnerg,PositionSameCompo(energyVal,:));
  
    
%% Plot all particles for a given energy (by step of 3) %% 
    PlotParticleTrajectory(electronsV0Neg,EnergyPartsIndices(energyVal):3:EnergyPartsIndices(energyVal+1)-1, 1:nrun)

%% Plot all particles with given initial components for fixed E %%
    disp('All particles with same initial velocity components for a given energy')
    disp(strcat('(V0R,V0Z) = (',num2str(PartInfoVneg(5,PositionSameCompo(energyVal,1))), ',' ...
                               ,num2str(PartInfoVneg(4,PositionSameCompo(energyVal,1))) , ')'))
    disp(strcat('E = ', num2str(E(energyVal)), ' eV'));
    PlotParticleTrajectory(electronsV0Neg,PositionSameCompo(energyVal,:),1:nrun/10)
    
%% Plot all components for a given particle position %%
    PlotParticleTrajectory(electronsV0Neg, PosAllCompoPerPart(energyVal,:), 1:nrun)
    

    
    
    
    
    
    
    
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
    
    
    
     
%% Compare Vn and purely normal V0 scan %%

    compoVal  = 6; % purely normal (r) initial velocity
    energyVal = 2; % energy value = E(energyVal)
    posVal    = 4; % initial electron position
    
    % Find all indices with same (vr,vz) for given value %
    for ii = 1:nPoints 
        
        PositionSameCompo(ii,:) = (1+nElectrons*(compoVal-1)) + (nElectrons*nComponents)*(ii-1) : (1+nElectrons*(compoVal-1) + nElectrons-1  ) + (nElectrons*nComponents)*(ii-1);
        
    end
    
    % find all indices for given energy value %
    for ii = 1:nPoints
        
        PositionSameEnergVn(ii,:) = (1+nElectrons*(ii-1)):(1+nElectrons*(ii-1)+nElectrons -1);
        
    end
    
    % find all components for same particle (psoition) %
    for ii =1:nPoints 
       
        for jj = 1:nComponents 
            
            PosAllCompoPerPart(ii,jj) = nElectrons*(jj-1)+posVal +nElectrons*nComponents*(ii-1); % eventually change nElec*nCompo = 72
            
        end 
    end
    

    PlotParticleTrajectory(electronsV0,PositionSameCompo(energyVal,:),1:nrun)     % electrons with V0 along r
    
    PlotParticleTrajectory(electronsVn, PositionSameEnergVn(energyVal,:), 1:nrun) % electrons with Vr (scan)
    
    

%% Plot energy ratio %% 


