    %% Process electron trajectories for electrons generated at electrodes %%

    path = '/scratch/sguincha';
    directory = '/scratch/sguincha/gt170_refurbished_4_6_25kV';
    %-----------------------------------------
    addpath(genpath(strcat(path,directory)));
    addpath(genpath('/home/sguincha/espic2d/matlab/'))
    addpath(genpath('/home/sguincha/SPE_Fall2022/Matlab/'))
    format long 

    electrons   = espic2dhdf5('resultfast.h5');
    
    %% Display fields %%
    dispespicFields(electrons)
    
    
    %% Extract extrude Geometry %%
    % Recall that the geometry is elliptical 
   
    