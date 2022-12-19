%% Add all paths %%
addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines

%% PROCESSING 
filename = 'H2Slanted_1e-12.h5';
ions     = espic2dhdf5(filename);

%%
ions.displaytotcurrevol_geom(length(ions.t2d)-500:10:length(ions.t2d))