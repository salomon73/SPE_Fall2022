%% UPDATE RESTART FILE %%

electrons_gauss_iiee_id = -1;     % no ion_induced for el (default)
electrons_gauss_material_id = 1; % stainless steel (default)
electrons_gauss_zero_vel = 1; % zero initial vel. (default)
electrons_gauss_neuttype_id = 1; % hydrogen 


ion_tracers_iiee_id    = 2;   % ions add el to species 2
ion_tracers_material_id = 1;  % ions impinging on 304 SS
ion_tracers_zero_vel = 1;  % zero initial vel. (default) 
ion_tracers_neuttype_id = 1;  % hydrogen 


electron_tracers_iiee_id = -1;     % el. don't add to any species
electron_tracers_material_id = 1;  % 304 SS (default)
electron_tracers_zero_vel = 0; % avg init. vel. 2eV
electron_tracers_neuttype_id = 1;  % hydrogen (default)

weight  = 5.8653e+04;
masses  = [9.10938300000000e-31, 9.10938300000000e-31, 1.67262192000e-27]';
weights = [5.8653e+04, 5.8653e+04, 5.8653e+04]';

Partfile_order = ['electrons_gauss.in', 'electron_tracers.in', 'ion_tracers.in'];

fid = H5F.open('restartfast.h5','H5F_ACC_RDWR','H5P_DEFAULT');
    h5disp('restartfast.h5')
    group1_id = H5G.open(fid, '/Parts');
    
        h5writeatt('restartfast.h5', '/Parts/', 'iiee_id', int32(electrons_gauss_iiee_id) )
        h5writeatt('restartfast.h5', '/Parts/', 'material_id', int32(electrons_gauss_material_id) )
        h5writeatt('restartfast.h5', '/Parts/', 'neuttype_id', int32(electrons_gauss_neuttype_id) )
        h5writeatt('restartfast.h5', '/Parts/', 'zero_vel', int32(electrons_gauss_zero_vel))

        h5writeatt('restartfast.h5', '/Parts/ 2', 'iiee_id', int32(electron_tracers_iiee_id) )
        h5writeatt('restartfast.h5', '/Parts/ 2', 'material_id', int32(electron_tracers_material_id) )
        h5writeatt('restartfast.h5', '/Parts/ 2', 'neuttype_id', int32(electron_tracers_neuttype_id) )
        h5writeatt('restartfast.h5', '/Parts/ 2', 'zero_vel', int32(electron_tracers_zero_vel) )
        
        h5writeatt('restartfast.h5', '/Parts/ 3', 'iiee_id', int32(ion_tracers_iiee_id) )
        h5writeatt('restartfast.h5', '/Parts/ 3', 'material_id', int32(ion_tracers_material_id) )
        h5writeatt('restartfast.h5', '/Parts/ 3', 'neuttype_id', int32(ion_tracers_neuttype_id) )
        h5writeatt('restartfast.h5', '/Parts/ 3', 'zero_vel', int32(ion_tracers_zero_vel) )
        h5writeatt('restartfast.h5', '/Parts/ 3', 'is_field', int32(0) )
    
    
        h5write('restartfast.h5', '/Parts/masses', masses)
        h5write('restartfast.h5', '/Parts/weights', weights)
        h5disp('restartfast.h5')


