% GT170 1MW geometry
% Contains the geometry of the first 1MW 170GHz CW TH1509 tube
% JPH : Aug. 4th, 2017
%
	z_shim              = 0.00;                      % Shim that is possibly used to lift the gyrotron

    z_cat_min_design    = 0.1335;                     % Coordinate of emitter at larger radius
    z_cat_design        = 0.136;                     % 1MW 170GHz first prototype (2017), cathode CENTER
	z_cav_design        = 0.4995;
    z_launcher_design   = 0.78848;                   % Launcher edge
    
    r_cat_min           = 0.05419;                   % largest emitter radius of first prototype 170GHz 1MW
    r_cat               = 0.0533693;                 % 1MW 170GHz first prototype (2017) AVERAGE radius, WARM dimension
    r_launcher          = 0.02068;                   % launcher radius

% Add shimming if necessary

    z_cat_min           = z_cat_min_design  + z_shim;
    z_cat               = z_cat_design      + z_shim;
	z_cav               = z_cav_design      + z_shim;
    z_launcher          = z_launcher_design + z_shim;

