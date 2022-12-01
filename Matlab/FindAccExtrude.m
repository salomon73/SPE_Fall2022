%% Find potential bias for ions in extrude geometry %%
%  to estimate the min required number of time steps %

format long

rgrid = ions.rgrid;
zgrid = ions.zgrid;

rpos_dn = 0.00880;
rpos_up = 0.01173;

zpos_lf = 0.3570;
zpos_rg = 0.4025;

z_bot = 0.3751;
r_bot = 0.007533;

r_top = 0.0122;
z_top = 0.3784;

diffr_top = abs(r_top-ions.rgrid);
diffz_top = abs(z_top-ions.zgrid);

diffr_bot = abs(r_bot-ions.rgrid);
diffz_bot = abs(z_bot-ions.zgrid);

diffr_dn = abs(rpos_dn-ions.rgrid);
diffr_up = abs(rpos_up-ions.rgrid);

diffz_lf = abs(zpos_lf-ions.zgrid);
diffz_rg = abs(zpos_rg-ions.zgrid);

[val_dn, Id_dn] = min(diffr_dn);
[val_up, Id_up] = min(diffr_up);


[val_lf, Id_lf] = min(diffz_lf);
[val_rg, Id_rg] = min(diffz_rg);

[val_bot_r, Id_bot_r] = min(diffr_bot);
[val_bot_z, Id_bot_z] = min(diffz_bot);

[val_top_r, Id_top_r] = min(diffr_top);
[val_top_z, Id_top_z] = min(diffz_top);


pot = ions.pot;
pot = pot(:,:,end);

pot_up = pot(Id_top_r, Id_bot_z);
pot_dn = pot(Id_bot_r, Id_bot_z);


bias = abs(pot_up - pot_dn);
bias_up_down = abs(pot(length(rgrid), Id_bot_z) - pot(Id_bot_r, Id_bot_z) );

rb = 0.02355;
ra = 0.007533;

r_cloud = 0.009867;
r_bot_cloud = 0.008333;
diff_r_cloud = abs(r_cloud - rgrid);
[val_r_cloud, Id_r_cloud] = min(diff_r_cloud); 
diff_rbot_cloud = abs(r_bot_cloud - rgrid);
[val_rbot_cloud, Id_rbot_cloud] = min(diff_rbot_cloud); 


bias_cloud = abs(pot(Id_r_cloud, Id_bot_z) - pot(Id_bot_r, Id_bot_z) );
bias_min   = abs(pot(Id_r_cloud, Id_bot_z) - pot(Id_rbot_cloud, Id_bot_z) );

mion = 3.347e-27;
Ekin_max  = bias_cloud * log(rgrid(Id_top_r)/ra)/ log(rb/ra)
vmax = sqrt(2*Ekin_max/mion*1.602e-19)
Ekin_min  = bias_min * log(rgrid(Id_top_r)/ra)/ log(rb/ra)
vmin = sqrt(2*Ekin_min/mion*1.602e-19)



% find number of time steps for simulation 
vbar = (vmin + vmax)/2;
dt   = 5e-12; 
DeltaT = r_top/vbar
Nrun   = DeltaT/dt
Nrun   =  round(Nrun*2)