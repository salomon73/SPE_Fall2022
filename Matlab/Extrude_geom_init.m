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
    
    r0 = 0.01;     %geom.r_0;
    z0 = geom.z_0;
    dr = 1e-8;
    
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

    % tangent vector 
    tgt1 = -zr*sin(tel);
    tgt2 =  rr*cos(tel);
    
    tgt = [tgt1; tgt2];
    
    % normal vector
    nor1 = -tgt2; % r component
    nor2 =  tgt1; % z component
   
    nor = [nor1; nor2]; 


    geomstruct = geom;
    figure
        subplot(1,2,1)
            contour(geomstruct.zgrid,geomstruct.rgrid,geomstruct.geomweight(:,:,1),[0 0],'k-', 'linewidth', 2)
            [~,rgrid]=meshgrid(geomstruct.zgrid,geomstruct.rgrid);
            rlims=rgrid(geomstruct.geomweight(:,:,1)<0);
            ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on        
            plot(zel, rel , 'r*', 'linewidth', 2)
            hold on
            quiver(zel,rel,...
                   nor2, nor1,...
                   'b');
            legend('electrode', '$e^-$' ,...
                        'Location','northwest','Interpreter','latex');
        subplot(1,2,2)
            hold on 
            plot(1e3*zcat, 1e3*rcat , 'k-', 'linewidth', 2)
            hold on 
            plot(1e3*zel, 1e3*rel , 'r*', 'linewidth', 2)
            ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            legend('electrode', '$e^-$','Location','northwest','Interpreter','latex');


 
                    
    %% GENERATE INPUT FOR EXTRUDE GEOMETRY %% 
   
    % Physical constants 
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    format long;
    
    % Energy characteristics 
    lowerBound = 0.1;
    upperBound = 10;
    nPoints    = 20;
    nEnergy    = nPoints;
    E = linspace(lowerBound, upperBound, nPoints);

    % Input electrons characteristics 
    nElectrons = length(tel);
    nblocks = nEnergy * nElectrons; 
    nparts  = nEnergy * nElectrons;
    
    
    % Velocity characteristics (components initialisation)
    VNorm = zeros(length(E), length(nor(:,1)), length(nor(1,:))); 
    
    for kk = 1:length(E)
        
        VNorm(kk,:,:) = sqrt(2*E(kk)*e/m) * nor(:,:);
        
    end

    
    % File characteristics 
    
    TextFileVn    = 'V0Extrude';
    
    
    % File initialisation (Vn)
        
    fileId = fopen(strcat(TextFileVn,'.txt'),'w');
    fprintf(fileId,'%s %d \n', 'nblocks =' , nblocks );
    fprintf(fileId,'%s %d \n', 'nparts ='  , nparts );
    fprintf(fileId,'//parts\n');
    for ii =1:length(E)
        for jj =1:length(tel)
            fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', rel(jj)+dr, 0.0, zel(jj), VNorm(ii,1,jj), 0.0 ,VNorm(ii,2,jj));
        end
    end
    fclose(fileId);
    
    % Write particle informations in variable file
    PartInfoVn   = zeros(6,nparts);
    compteur = 0;
        
    for ii =1:length(E)

        for jj =1: nElectrons
            
                compteur = compteur + 1; 
                PartInfoVn(1,compteur) = compteur;  
                PartInfoVn(2,compteur) = rel(jj)+dr;
                PartInfoVn(3,compteur) = zel(jj);
                PartInfoVn(4,compteur) = VNorm(ii,1,jj);
                PartInfoVn(5,compteur) = VNorm(ii,2,jj);
                PartInfoVn(6,compteur) = E(ii);

        end
    end
    save(strcat('PartInfos_',num2str(nEnergy),'_',num2str(nElectrons),'.mat'),'PartInfoVn');
