%% Processing the data for self sustained iiee in several TREX geometries %%

addpath(genpath('/home/sguincha/espic2d/matlab'))
ions = espic2dhdf5('H2Slanted_1e-12.h5') 


%% display current collected in different regions
ions.displaytotcurrevol(1:10:length(ions.t2d))

%% display surface current density collected in different regions
ions.displaySurfFlux(length(ions.t2d))

