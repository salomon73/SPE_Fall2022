function AddAllPaths()
path = '/scratch/sguincha';
directory = '/SPE_HRES_RUN/ResultsScan3/'; % Change result folder name accordingly
directory2 = 'SPE_Fall2022/RunsElectronsSlices/Results/';
directory3 = '/scratch/sguincha/SPE_Fall2022/Results/';
directory4 = '/scratch/sguincha/gt170_refurbished_4_6_25kV';
%-----------------------------------------
addpath(genpath(strcat(path,directory)));
addpath(genpath(strcat(path,directory2)));
addpath(genpath(strcat(path,directory3)));
addpath(genpath('/home/sguincha/espic2d/matlab/'))
addpath(genpath('/home/sguincha/SPE_Fall2022/Matlab/'))
addpath(genpath(directory4))
format long 
end