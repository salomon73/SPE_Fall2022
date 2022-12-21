%% Processing the data for self sustained iiee in several TREX geometries %%

addpath(genpath('/home/sguincha/espic2d/matlab'))
ions = espic2dhdf5('H2Slanted_1e-12.h5') 


%% 

ions.displaytotcurrevol(1:10:length(ions.t2d))