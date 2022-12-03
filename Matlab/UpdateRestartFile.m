%% UPDATE RESTART FILE %%

electrons_gauss_iiee_d = -1;     % no ion_induced for el (default)
electrons_gauss_material_id = 1; % stainless steel (default)
electrons_gauss_zero_vel = true; % zero initial vel. (default)
electrons_gauss_neuttype_id = 1; % hydrogen 


ion_tracers_iiee_id    = 3;   % ions add el to species 3
ion_tracers_material_id = 1;  % ions impinging on 304 SS
ion_tracers_zero_vel = true;  % zero initial vel. (default) 
ion_tracers_neuttype_id = 1;  % hydrogen 


electron_tracers_iiee_id = -1;     % el. don't add to any species
electron_tracers_material_id = 1;  % 304 SS (default)
electron_tracers_zero_vel = false; % avg init. vel. 2eV
electron_tracers_neuttype_id = 1;  % hydrogen (default)

Partfile_order = ['electrons_gauss.in', 'electron_tracers.in', 'ion_tracers.in'];

fid = H5F.open('restartfast.h5','H5F_ACC_RDWR','H5P_DEFAULT');
    h5disp('restartfast.h5')
    group1_id = H5G.open(fid, '/Parts');
        h5write('/Parts/iiee_id', electron_tracers_iiee_id);
        h5write('/Parts/material_id', electron_tracers_material_id);
        h5write('/Parts/neuttype_id', electron_tracers_neuttype_id);
        h5write('/Parts/zero_vel', electron_tracers_zero_vel);
        
        
        
%         h5write('/Parts/ 2/iiee_id', electrons_gauss_iiee_d);
%         h5write('/Parts/ 2/material_id', electrons_gauss_material_id);
%         h5write('/Parts/ 2/neuttype_id', electrons_gauss_neuttype_id);
%         h5write('/Parts/ 2/zero_vel', electrons_gauss_zero_vel);
%                         
%         
%         h5write('/Parts/ 3/iiee_id', ion_tracers_iiee_id);
%         h5write('/Parts/ 3/material_id', ion_tracers_material_id);
%         h5write('/Parts/ 3/neuttype_id', ion_tracers_neuttype_id);
%         h5write('/Parts/ 3/zero_vel', ion_tracers_zero_vel);
%                         
% 
