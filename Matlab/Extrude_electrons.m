    %% Process electron trajectories for electrons generated at electrodes %%

    path = '/scratch/sguincha';
    directory = '/scratch/sguincha/gt170_refurbished_4_6_25kV';
    %-----------------------------------------
    addpath(genpath(strcat(path,directory)));
    addpath(genpath('/home/sguincha/espic2d/matlab/'))
    addpath(genpath('/home/sguincha/SPE_Fall2022/Matlab/'))
    format long 

    electrons   = espic2dhdf5('resultfast.h5');
    geom        = electrons;
    
    %% Display fields %%
    dispespicFields(electrons)
    
    
    %% Extract extrude Geometry %%
    % Draw the cathode
    % ellipse equation (z-z0)^2/a^2 + (r-r0)^2/b^2 = 1 
    
    r0 = 0.01; %geom.r_0;
    z0 = geom.z_0;
    
    zr = geom.z_r; % a semi major axis
    rr = geom.r_r; % b semi minor axis
    nppts = 1000;  % # of points to get a smooth ellipse
    tcat  = pi/2:0.01:3*pi/2;
    rcat  = r0 + rr*cos(tcat);
    zcat  = z0 + zr*sin(tcat);  
 
    

    
    % initialise nElectrons on the internal elliptical surface
    nElectrons = 20;
    tel  = pi/2:pi/nElectrons:3*pi/2;
    rel  = r0 + rr*cos(tel);
    zel  = z0 + zr*sin(tel); 
   
    figure
        hold on 
        plot(1e3*zcat, 1e3*rcat , 'k-', 'linewidth', 2)
        hold on 
        plot(1e3*zel, 1e3*rel , 'r*', 'linewidth', 2)
        ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
        xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
        set (gca, 'fontsize', 20)
        legend('electrode', '$e^-$','Location','northwest','Interpreter','latex');
        


    geomstruct = geom;
    figure
            contour(geomstruct.zgrid,geomstruct.rgrid,geomstruct.geomweight(:,:,1),[0 0],'k-', 'linewidth', 2)
            [~,rgrid]=meshgrid(geomstruct.zgrid,geomstruct.rgrid);
            rlims=rgrid(geomstruct.geomweight(:,:,1)<0);
            ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on        
            plot(zel, rel , 'r*', 'linewidth', 2)
%             plot(1e3*Points(2,1:end),1e3*Points(1,1:end),'r*')
%             hold on 
%             quiver(1e3*Points(2,2:3:end),1e3*Points(1,2:3:end),...
%                    scalefactor*V0(1,2:3:end,end,1), scalefactor*V0(1,2:3:end,end,2),...
%                    'b', 'Autoscale', 'off');
            legend('electrode', '$e^-$' ,...
                        'Location','northwest','Interpreter','latex');

 
                    
    %% GENERATE INPUT FOR EXTRUDE GEOMETRY %% 
        