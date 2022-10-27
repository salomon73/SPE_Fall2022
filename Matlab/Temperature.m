function velocity  = Temperature()

%% Define the constants %%
kB = (25.7/298)*1e-3;      % eV/K
m  = 9.10938300000000e-31; % electron mass
e  = 1.60217662000000e-19; % J/eV

%% Energy range for T range %%
velocity.E  = linspace(0.1, 15, 10);
velocity.T  = 1/kB*velocity.E;
velocity.v  = sqrt((2*e/m)*velocity.E);
end