% This script tests the coefficients of the matrix found by
% 'B_Mat88_Cryogenic_170
%
% JHE, July 31st, 2017
%

    I= 100*[0.57712  , 0.68596 ,  1.07780  , 1.07780] ; % Nominal currents used on MAy 29th, 2017 to apply first SP pulses
    

%     r_cat_min			= 0.05419;                  % largest emitter radius of first prototype 170GHz 1MW
%     r_launcher          = 0.02068;                  % launcher radius
% 
%     z_cat_min 	        = 0.1335;    
%     z_launcher          = 0.78848;


%  Include GT170 Geometry
    
    GT170_1MW_Geometry
    
%
% Taking parameters defined with on-axis fields
    %gamma = [0.002933681588584, -0.006273546452915,  -0.007626766134967, 0.009858220586519];

% Taking more stringent condition (on axis field at gun, off-axis field at launcher
%   gamma = [ 0.00293370270241, -0.00627361773021,   -0.00762699351991 , 0.00987214759867];

% Taking parameters defined with potential vector (defined as delta in
% B_Matrix_Cryogenic_170, with another factor 1000 omitted
    gamma = [0.081474637274171, -0.174106309671082, -0.203575059328030,  0.260894824884893];

    
    z = linspace(z_launcher-0.01,z_launcher+0.01,41);
    r = linspace(0.015, 0.030,16);
    [Z,R] = meshgrid(z,r);
    
    [aphi_cat] = B_Ellip_Cryogenic_170('aphi','cryogenic',I,r_cat_min,z_cat_min);
    raphi_cat  = r_cat_min * aphi_cat;
    
    [Aphi]     = B_Ellip_Cryogenic_170('aphi','cryogenic',I,R,Z);
    RAphi      = R.*Aphi;
    
    figure(1) 
    contour(Z,R,RAphi,[raphi_cat raphi_cat],'--','LineWidth',2)
    
    
    hold on
    
    I(4) = I(4) + 10;
    [aphi_cat] = B_Ellip_Cryogenic_170('aphi','cryogenic',I,r_cat_min,z_cat_min);
    raphi_cat  = r_cat_min * aphi_cat;
    
    [Aphi]     = B_Ellip_Cryogenic_170('aphi','cryogenic',I,R,Z);
    RAphi      = R.*Aphi;
    
    figure(1) 
    contour(Z,R,RAphi,[raphi_cat raphi_cat],'-.','LineWidth',2)
    
    
    I4lim = - (gamma(1) * I(1) +gamma(2)*I(2) +gamma(3)*I(3)) / gamma(4)
    
    I(4) = I4lim;
    [aphi_cat] = B_Ellip_Cryogenic_170('aphi','cryogenic',I,r_cat_min,z_cat_min);
    raphi_cat  = r_cat_min * aphi_cat;
    
    [Aphi]     = B_Ellip_Cryogenic_170('aphi','cryogenic',I,R,Z);
    RAphi      = R.*Aphi;
    
    figure(1) 
    contour(Z,R,RAphi,[raphi_cat raphi_cat],'LineWidth',2)
    plot(z_launcher,r_launcher,'r*')
   