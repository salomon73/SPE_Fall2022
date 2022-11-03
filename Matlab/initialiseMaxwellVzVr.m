    %% Add all paths containing outputs %%
    AddAllPaths();
    
    %% Define the constants %%
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    
    %% GT170 GYROTRON GEOMETRY %%
    geom = espic2dhdf5('result_46dA_25kv_10mubar.h5');
    
    %% Create geom structure %% 

    lowerBound = 0.1;
    upperBound = 20;
    nPoints = 20;
    cathode = geom.spl_bound.boundary(1);
    LowPts  = 15;
    UpPts   = 65;
    nblocks = nPoints*(UpPts-LowPts);
    nparts  = nblocks;
    dr      = 1e-6;
    nElectrons = UpPts-LowPts;
    
    f   = cathode.fun;
    E   = linspace(lowerBound, upperBound, nPoints);
    Points = zeros(UpPts-LowPts,2);    

    Tang   = zeros(UpPts-LowPts,2);
    Norm   = zeros(UpPts-LowPts,2);
    VNorm  = zeros(length(E), UpPts-LowPts,2);
    df = fnder(f,1);
    TextFileVnorm = 'VNorm';
    TextFileV0    = 'V0';
    
    for ii = 0:UpPts-LowPts-1
        
       pts = spval(f,f.knots(ii+LowPts));
       Points(ii+1,1) = pts(1);
       Points(ii+1,2) = pts(2);
       tgt = spval(df,df.knots(ii+LowPts-1));
       Tang(ii+1,1) =  tgt(1);
       Tang(ii+1,2) =  tgt(2);
       Norm(ii+1,1) = -tgt(2);
       Norm(ii+1,2) =  tgt(1);
    end

    for jj =1:length(Norm(:,1))

        Norm(jj,:) = 1/norm(Norm(jj,:),2)*Norm(jj,:);   
    end
    
    for kk = 1:length(E)
        
        VNorm(kk,:,:) = sqrt(2*E(kk)*e/m) * Norm(:,:);
        
    end

    LowInd  = LowPts-1;
    HighInd = UpPts-1;
    
    figure
        hold on 
        plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
        ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
        xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
        set (gca, 'fontsize', 20)
        hold on 
        plot(Points(1:3:end,1),Points(1:3:end,2),'r*')
        hold on 
        quiver(Points(1:3:end,1),Points(1:3:end,2), Norm(1:3:end,1), Norm(1:3:end,2), 'b')
        axis equal
        
    
    fileId = fopen(strcat(TextFileVnorm,'.txt'),'w');
    fprintf(fileId,'%s %d \n', 'nblocks =' , nblocks );
    fprintf(fileId,'%s %d \n', 'nparts =' , nparts );
    fprintf(fileId,'//parts\n');
    for ii =1:length(E)
        for jj =1:UpPts-LowPts
            fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', Points(jj,2)+dr, 0.0, Points(jj,1), VNorm(ii,jj,2), 0.0 ,VNorm(ii,jj,1));
        end
    end
    fclose(fileId);
    
    %%  Plot normal vectors times velocity for given energy values %%
        figure
        subplot(2,2,1)
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
            hold on 
            quiver(Points(1:3:end,1),Points(1:3:end,2), VNorm(1,1:3:end,1)', VNorm(1,1:3:end,2)','b')
            legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(1)), ' eV'), 'Location','northwest','Interpreter','latex');
            axis equal
        subplot(2,2,2)
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
            hold on 
            quiver(Points(1:3:end,1),Points(1:3:end,2), VNorm(5,1:3:end,1)', VNorm(5,1:3:end,2)', 'b')
            legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(5)), ' eV'), 'Location','northwest','Interpreter','latex');
            axis equal
        subplot(2,2,3)
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
            hold on 
            quiver(Points(1:3:end,1),Points(1:3:end,2), VNorm(10,1:3:end,1)', VNorm(10,1:3:end,2)', 'b')
            legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(10)), ' eV'), 'Location','northwest','Interpreter','latex');
            axis equal
        subplot(2,2,4)
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
            quiver(Points(1:3:end,1),Points(1:3:end,2), VNorm(end,1:3:end,1)', VNorm(end,1:3:end,2)', 'b')
            legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(end)), ' eV'), 'Location','northwest','Interpreter','latex');
            axis equal
            
    %% Plot normal velocity vectors for min and max energies %%
    scaleFactor = sqrt(E(end)/E(1));
    scale       = 2e-8;
         figure
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:5:end,1),Points(1:5:end,2),'r.', 'markersize',12)
            hold on 
            quiver(Points(1:5:end,1),Points(1:5:end,2), VNorm(1,1:5:end,1)', VNorm(1,1:5:end,2)','b')
            hold on 
            plot(Points(2:5:end,1),Points(2:5:end,2),'g.', 'markersize',12)
            hold on 
            quiver(Points(2:5:end,1),Points(2:5:end,2), scale*VNorm(end,2:5:end,1)', scale*VNorm(end,2:5:end,2)','r','Autoscale', 'off')
            legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(1)), ' eV'), '$e^-$', strcat('$\mathbf{v_0}$',' : E =',num2str(E(end)), ' eV'), 'Location','northwest','Interpreter','latex');
            axis equal
    
    
    %% Influence of normal and tangential velocity components %% 

    format long
    nComponents = 7; % number of values for velocity
    nppts = UpPts-LowPts;
    vR    = zeros(length(E),nComponents); 
    vZ    = zeros(length(E),nComponents); 

    % vR and vZ for constant E 
    for kk=1:length(E)    
        vRmax = sqrt(2*E(kk)/m*e);
        vR(kk,:)    = linspace(0,vRmax, nComponents);

        for j = 1 :length(vR(kk,:))
            if 2*E(kk)/m*e - vR(kk,j)^2 < 1e-3
                vZ(kk,j) = 0.0;
            else 
                vZ(kk,j) = sqrt(2*E(kk)/m*e - vR(kk,j)^2);
            end
        end
    end    
    
    V0 = zeros(length(E),nppts,nComponents,2); % nEnergy * nParticles * nVelocityValues * 2 (r,z) components for each vector
    
    for ii = 1: length(E)
       
        for jj =1:nComponents
           V0(ii,:,jj,1) = vZ(ii,jj);
           V0(ii,:,jj,2) = vR(ii,jj);
        end
        
    end

    nblockV0 = length(E) * nComponents * nppts;
    npartsV0   = nblockV0;
    fileId = fopen(strcat(TextFileV0,'.txt'),'w');
    fprintf(fileId,'%s %d \n', 'nblocks =' , nblockV0 );
    fprintf(fileId,'%s %d \n', 'nparts =' , npartsV0 );
    fprintf(fileId,'//parts\n');
    
    for ii =1 : length(E)
        for jj =1 : nComponents
            for kk =1 : nppts
                fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', Points(kk,2)+dr, 0.0, Points(kk,1), V0(ii,kk,jj,1), 0.0 ,V0(ii,kk,jj,2));
            end
        end
    end
    
    fclose(fileId);
    
  

        figure
            hold on 
            plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
            ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on 
            plot(Points(1:3:end,1),Points(1:3:end,2),'r*')
            hold on 
            quiver(Points(1:3:end,1),Points(1:3:end,2), V0(1,1:3:end,end,1)', V0(1,1:3:end,end,2)', 'b')
            axis equal

            
    %% Enumerate and label the particles for trajectories processing %%

    PartInfoV0 = zeros(6,npartsV0);
    PartInfoVn = zeros(6,nparts);
    compteur   = 0;
    
    % All informations about particles initialised with V0 % 
    for ii = 1:length(E)

       for kk = 1: nComponents

           for jj = 1: nElectrons
               
                compteur = compteur + 1; 
                PartInfoV0(1,compteur) = compteur;  
                PartInfoV0(2,compteur) = Points(jj,2);
                PartInfoV0(3,compteur) = Points(jj,1);
                PartInfoV0(4,compteur) = V0(ii,jj,kk,1);
                PartInfoV0(5,compteur) = V0(ii,jj,kk,2);
                PartInfoV0(6,compteur) = E(ii);
           end

       end

    end

    compteur = 0;
    % All informations about particles initialised with VN % 
    for ii =1:length(E)

        for jj =1: nElectrons
            
                compteur = compteur + 1; 
                PartInfoVn(1,compteur) = compteur;  
                PartInfoVn(2,compteur) = Points(jj,2);
                PartInfoVn(3,compteur) = Points(jj,1);
                PartInfoVn(4,compteur) = VNorm(ii,jj,1);
                PartInfoVn(5,compteur) = VNorm(ii,jj,2);
                PartInfoVn(6,compteur) = E(ii);

        end
    end

   save('PartInfosGT170.mat','PartInfoV0','PartInfoVn');
%==============================================================================================================
%==============================================================================================================
    
    %% TREX GEOMETRY %%
    geomTRex = espic2dhdf5('resultrestart_5e-12.h5');
    %dispespicFields(geomTRex)
    
    %% Define the constants %%
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    
 
    %% Energy scan %% 
    format long 
    
    % Energy parameters % 
    lowerBound = 0.1;
    upperBound = 20;
    nPoints = 10; % # values to scan E 
    E = linspace(lowerBound, upperBound, nPoints);
    
    % geom parameters %
    zA = 0.25;
    zB = 0.45;
    rA = geomTRex.r_a;
    rB = geomTRex.r_b;
    dr = 1e-6; % [m] to particle existence
    Zlim   = [0.30 0.43];
    
    % particles parameters %
    nElectrons = 12;
    PointsZ = linspace(Zlim(1), Zlim(2), nElectrons);
    PointsR = rA+dr;
    Points  = [ PointsR * ones(size(PointsZ));PointsZ ];
    
    % Numerical parameters %
    nblocks = nPoints * nElectrons;
    nparts  = nPoints * nElectrons;
    TRexTxtVNorm = 'vNormTRex';
    TRexTxtV0    = 'v0TRex';
    
    
    Norm  = zeros(nPoints, nElectrons, 2);
    VNorm = zeros(nPoints, nElectrons, 2);
    for ii =1:nPoints 
        for jj =1:nElectrons
            Norm(ii,jj,:) = [0 1]; %vecteur normal (0 1) 
        end
        VNorm(ii,:,:) = sqrt(2*E(ii)*e/m) * Norm(ii,:,:);
    end
    
   
    % Normal direction velocity scan for different E %
    
    fileId = fopen(strcat(TRexTxtVNorm,'.txt'),'w');
    fprintf(fileId,'%s %d \n', 'nblocks =' , nblocks );
    fprintf(fileId,'%s %d \n', 'nparts ='  , nparts );
    fprintf(fileId,'//parts\n');
    for ii =1:nPoints
        for jj =1:nElectrons
            fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', Points(1,jj), 0.0, Points(2,jj), VNorm(ii,jj,2), 0.0 ,VNorm(ii,jj,1));
        end
    end
    fclose(fileId);
    
    
    
    % Influence of velocity directions & 
    
    nComponents = 6; % number of values for velocity
    vR    = zeros(length(E),nComponents); 
    vZ    = zeros(length(E),nComponents); 
    
    % vR and vZ for constant E 
    for kk=1:length(E)    
        vRmax = sqrt(2*E(kk)/m*e);
        vR(kk,:)    = linspace(0,vRmax, nComponents);

        for j = 1 :length(vR(kk,:))
            if 2*E(kk)/m*e - vR(kk,j)^2 < 1e-3
                vZ(kk,j) = 0.0;
            else 
                vZ(kk,j) = sqrt(2*E(kk)/m*e - vR(kk,j)^2);
            end
        end
    end    
    

    V0 = zeros(length(E),nElectrons,nComponents,2); % nEnergy * nParticles * nVelocityValues * 2 (r,z) components for each vector
    
    for ii = 1: length(E)
       
        for jj =1:nComponents
           V0(ii,:,jj,1) = vZ(ii,jj);
           V0(ii,:,jj,2) = vR(ii,jj);
        end
        
    end
    
    nblockV0 = length(E) * nComponents * nElectrons;
    npartsV0   = nblockV0;
    fileId = fopen(strcat(TRexTxtV0,'.txt'),'w');
    fprintf(fileId,'%s %d \n', 'nblocks =' , nblockV0 );
    fprintf(fileId,'%s %d \n', 'nparts =' , npartsV0 );
    fprintf(fileId,'//parts\n');
    
    for ii =1 : length(E)
        for jj = 1:nComponents
            for kk = 1:nElectrons
                fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', Points(1,kk), 0.0, Points(2,kk), V0(ii,kk,jj,1), 0.0 ,V0(ii,kk,jj,2));
            end
        end
    end
    
    fclose(fileId);
    
    %% Plot initial condition for TREX scan %% 
    t=2*pi*linspace(0,1,100);
    
    geomstruct = geomTRex;
    scalefactor = 3e-5;
    figure
        subplot(1,2,1)
            contour(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3,geomstruct.geomweight(:,:,1),[0 0],'k-', 'linewidth', 2)
            [~,rgrid]=meshgrid(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3);
            rlims=rgrid(geomstruct.geomweight(:,:,1)<0);
            ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on
            plot(1e3*Points(2,1:end),1e3*Points(1,1:end),'r*')
            hold on 
            quiver(1e3*Points(2,2:3:end),1e3*Points(1,2:3:end),...
                   scalefactor*V0(1,2:3:end,end,1), scalefactor*V0(1,2:3:end,end,2),...
                   'b', 'Autoscale', 'off');
            legend('electrode', '$e^-$','$\mathbf{v_0} = v_0\mathbf{e_r}$',...
                        'Location','northwest','Interpreter','latex');


        subplot(1,2,2)
            contour(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3,geomstruct.geomweight(:,:,1),[0 0],'k-', 'linewidth', 2)
            [~,rgrid]=meshgrid(geomstruct.zgrid*1e3,geomstruct.rgrid*1e3);
            rlims=rgrid(geomstruct.geomweight(:,:,1)<0);
            ylabel('$R$ [mm]', 'interpreter', 'latex','Fontsize', 22)
            xlabel('$Z$ [mm]', 'interpreter', 'latex', 'Fontsize', 22)
            set (gca, 'fontsize', 20)
            hold on
            plot(1e3*Points(2,1:end),1e3*Points(1,1:end),'r*')
            for ii = 1:nComponents 
            hold on
            quiver(1e3*Points(2,2:3:end),1e3*Points(1,2:3:end), scalefactor*V0(1,2:3:end,ii,1),...
                   scalefactor*V0(1,2:3:end,ii,2),...
                   'b', 'Autoscale', 'off');
            end
                 legend('electrode', '$e^-$','$\mathbf{v_0} = v_r\mathbf{e_r}+v_z\mathbf{e_z}$',...
                        'Location','northwest','Interpreter','latex');

        
%% Enumerate and label the particles for trajectories processing %%

    PartInfoV0 = zeros(6,npartsV0);
    PartInfoVn = zeros(6,nparts);
    compteur   = 0;
    
    % All informations about particles initialised with V0 % 
    for ii = 1:length(E)

       for kk = 1: nComponents

           for jj = 1: nElectrons
               
                compteur = compteur + 1; 
                PartInfoV0(1,compteur) = compteur;  
                PartInfoV0(2,compteur) = Points(1,jj);
                PartInfoV0(3,compteur) = Points(2,jj);
                PartInfoV0(4,compteur) = V0(ii,jj,kk,1); % vZ
                PartInfoV0(5,compteur) = V0(ii,jj,kk,2); % vR
                PartInfoV0(6,compteur) = E(ii);
           end

       end

    end

    compteur = 0;
    % All informations about particles initialised with VN % 
    for ii =1:length(E)

        for jj =1: nElectrons
            
                compteur = compteur + 1; 
                PartInfoVn(1,compteur) = compteur;  
                PartInfoVn(2,compteur) = Points(1,jj);
                PartInfoVn(3,compteur) = Points(2,jj);
                PartInfoVn(4,compteur) = VNorm(ii,jj,1);
                PartInfoVn(5,compteur) = VNorm(ii,jj,2);
                PartInfoVn(6,compteur) = E(ii);

        end
    end
    
   save('PartInfos.mat','PartInfoV0','PartInfoVn');
